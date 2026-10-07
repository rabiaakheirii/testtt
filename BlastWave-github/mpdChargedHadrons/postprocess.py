#!/usr/bin/env python3

##
# @file postprocess.py
# @author Viktar Kireyeu
# @version 1.1
# @date    2024
# @copyright Do not modify a single character without the author's permission.
#
# @brief  Postprocessing macro for the MpdRoot 'nuclei' wagon.
#
# If the PTCORR file is set, then instead of the post-analysis the new
# ROOT-file with the \f$ p_{T} \f$ corrections 1d and 2d profiles will be produced.
#
# Usage: postprocess [-h] [-i INPUT] [-e EFFICIENCIES] [-s SETTINGS] [-o OUTPUT]
#                    [-d DIR] [-r REPORT] [-p PTCORRECTIONS] [-t TABLE] [--qa]
#                    [--png] [--dedx] [--skip] [--clear]
# 
# MpdNuclei wagon data processing program
# 
# Options:
# <pre>
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Input ROOT file
#   -e EFFICIENCIES, --eff EFFICIENCIES
#                         Efficiencies ROOT file
#   -s SETTINGS, --settings SETTINGS
#                         Settings JSON file
#   -o OUTPUT, --output OUTPUT
#                         Output spectra ROOT-file
#   -d DIR, --dir DIR     Output directory for histograms
#   -r REPORT, --report REPORT
#                         Output report file
#   -p PTCORRECTIONS, --ptcorr PTCORRECTIONS
#                         Make file for pt corrections
#   -t TABLE, --table TABLE
#                         Output file with dndy info for fits etc (JSON)
#   --qa                  Produce extended QA plots
#   --png                 Produce PNG-plots onstead of PDF
#   --dedx                Produce histograms for the case PID = dE/dx only
#   --skip                Skip already triggered functions
#   --clear               Clear the output directory before populating
# 
# If the EFFICIENCIES file is not set, then the INPUT will be used instead.
# </pre>
# 
# The new subdirectories will be created in the output directory (DIR):
# <pre>
# ├── contamination
# │   ├── pid
# │   ├── secondaries
# │   └── tof
# ├── dca
# ├── efficiency
# │   ├── pid
# │   ├── tof
# │   └── tpc
# └── results
#     ├── coal
#     ├── dndy
#     │   ├── ap
#     │   ├── d
#     │   ├── He3
#     │   ├── He4
#     │   ├── km
#     │   ├── kp
#     │   ├── p
#     │   ├── pim
#     │   ├── pip
#     │   └── t
#     ├── phasespace
#     │   ├── corrected
#     │   └── uncorrected
#     └── pt
#         ├── ap
#         ├── d
#         ├── He3
#         ├── He4
#         ├── km
#         ├── kp
#         ├── p
#         ├── pim
#         ├── pip
#         └── t
#</pre>
# 
# Changelog
# Version 1.1:
# - Bug fixes
# - New subroutines for the different PID methods
# - Multiprocessing for the pT-spectra fits
# - Documentation update
# 
# Version 1.0:
# - Initial version


import sys
try:
    import os
    import json
    import math
    import argparse
    import fileinput
    from array import array
    from pathlib import Path
    import random

    import ROOT
    from ROOT import TCanvas, TFile, TLatex, TLegend, TPad, TRatioPlot, TGraph, TGraphErrors, TF1, TMath, TObject, TPaveText
    from ROOT import gROOT, gSystem, gPad

    from pylatex import Document, Section, Subsection, Math, NoEscape, Package, Command
    import pylatex.config as cf

    import multiprocessing

except ModuleNotFoundError as err:
    sys.exit(err)

parser = argparse.ArgumentParser(
                    prog        = 'postprocess',
                    description = 'MpdNuclei wagon data processing program',
                    epilog      = 'If the EFFICIENCIES file is not set, then the INPUT will be used instead.')
parser.add_argument('-i', '--input',    metavar = 'INPUT',         help = 'Input ROOT file',                 default = 'input/train_output/paper_mpdpid10.root')
parser.add_argument('-e', '--eff',      metavar = 'EFFICIENCIES',  help = 'Efficiencies ROOT file')
parser.add_argument('-s', '--settings', metavar = 'SETTINGS',      help = 'Settings JSON file',              default = 'input/train_output/mpdpid_10prc.json')
parser.add_argument('-o', '--output',   metavar = 'OUTPUT',        help = 'Output spectra ROOT-file',        default = 'input/postprocess_mpdpid10.root')
parser.add_argument('-d', '--dir',      metavar = 'DIR',           help = 'Output directory for histograms', default = 'pyplots')
parser.add_argument('-r', '--report',   metavar = 'REPORT',        help = 'Output report file')
parser.add_argument('-p', '--ptcorr',   metavar = 'PTCORRECTIONS', help = 'Make file for pt corrections')
parser.add_argument('-t', '--table',    metavar = 'TABLE',         help = 'Output file with dndy info for fits etc (JSON)')
parser.add_argument('--qa',             help = 'Produce extended QA plots',                        action = 'store_true')
parser.add_argument('--png',            help = 'Produce PNG-plots onstead of PDF',                 action = 'store_true')
parser.add_argument('--dedx',           help = 'Produce histograms for the case PID = dE/dx only', action = 'store_true')
parser.add_argument('--skip',           help = 'Skip already triggered functions',                 action = 'store_true')
parser.add_argument('--clear',          help = 'Clear the output directory before populating',     action = 'store_true')


rapidity_bins = [[-1.0, -0.5], [-0.5, 0.5], [0.5, 1.0]] # Selected rapidity bin for the analysis
SMALL_DELTA = 0.00001
LW       =  5 # Line width
MS_MC    = 20 # Marker style: full circles
MS_UNCOR = 24 # Marker style: open circles
MS_COR   = 22 # Marker style: full triangles
LC_MC    =  1 # Line color: black
LC_UNCOR =  2 # Line color: red
LC_COR   =  3 # Line color: green

STATUS_FILE = 'postprocess.status'

# This is the main function, all other subroutines are called from here
def main():
    global EXT          # File extension for the produced plots: 'pdf', 'png', maybe something else
    global Particles    # Global list with the particles definitions from the settings file
    global Centrality   # Global list with the centralities definitions from the settings file
    global SYS          # Global variable with the colliding system title
    global DEDX         # Global switch for the case 'PID = dE/dx only': 
                        #   1 = make such plots; 0 = make plots for the combined PID case
    global CURRENT_ST   # Global list of the current functions calls counters
    global SKIP         # Global switch for skipping already triggered (during previous runs) functions
    global INP_FILE     # Input ROOT file
    global EFF_FILE     # Efficiencies ROOT file
    global OUT_FILE     # Output spectra ROOT-file
    global TBL_FILE     # Output file (JSON) to dump the BW fits and other information
    global PID_MODE     # The PID mode: 0 = evPID wagon, 2 = MpdPid class

    # ~ CURRENT_ST = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    CURRENT_ST = [0] * 20
    SYS = 'Bi+Bi, #sqrt{s_{NN}} = 9.2 GeV'

    # Command line arguments parsing
    args = parser.parse_args()
    input_file      = args.input           # Input file
    efficiency_file = args.eff             # Efficiencies file
    settings_file   = args.settings        # Settings file
    spectra_file    = args.output          # Output file for the spectra
    outdir          = args.dir             # Output directory to store plots
    extended        = args.qa              # Switch for the additional efficiencies plots
    report_file     = args.report          # Switch for the report file generation
    TBL_FILE        = args.table           # Output file for the BW fits info
    if args.eff is None:                   # If efficiencies file is not set
        efficiency_file = input_file       # we can use the input file

    EXT   = 'png' if args.png   else 'pdf' # pdf/png switch
    DEDX  = True  if args.dedx  else False # Switch for the 'PID = dE/dx only' case
    SKIP  = True  if args.skip  else False # Switch to skip previously triggered functions

    make_ptcorr = False
    if args.ptcorr is not None:            # If PTCOR file is set
        corrections_file = args.ptcorr     # then the program will only produce
        make_ptcorr = True                 # this file without furter analysis

    print(f'Input file: {input_file}')
    check_file_exist(input_file)                  # Check if the input file exists
    print(f'Settings file: {settings_file}')
    check_file_exist(settings_file)               # Check if the settings file exists

    if(make_ptcorr):
        print('*** pT corrections mode ***')
        print(f'pT corrections file: {corrections_file}')
    else:
        print(f'Efficiencies file: {efficiency_file}')
        check_file_exist(efficiency_file)         # Check if the efficiencies file exists
        print(f'Spectra file: {spectra_file}')
        if report_file is not None:
            print(f'Report: {report_file}')
        if TBL_FILE is not None:
            print(f'JSON table: {TBL_FILE}')
            if os.path.exists(TBL_FILE):
                os.system(f'rm -f {TBL_FILE}')    # Remove the JSON 'table' file if it exists
        if args.clear:                            # The 'dry' run case -- clear everything before proceed
            os.system(f'rm -rf {outdir}')         # Clear remnants from previous runs
            os.system(f'rm -vf {STATUS_FILE}')
            os.system(f'rm -vf {spectra_file}')
            os.system(f'mkdir -pv {outdir}/dca')
            os.system(f'mkdir -pv {outdir}/results/{{pt,dndy,coal}}')
            os.system(f'mkdir -pv {outdir}/results/phasespace/{{uncorrected,corrected}}')
        if extended:                              # The extended output case -- create the needed sub-directories
            os.system(f'mkdir -pv {outdir}/efficiency/{{tpc,tof,pid}}')
            os.system(f'mkdir -pv {outdir}/contamination/{{secondaries,pid,tof}}')
            if DEDX:                              # Same for the 'PID = dE/dx only' case
                os.system(f'mkdir -pv {outdir}/efficiency/pid_dedx')
                os.system(f'mkdir -pv {outdir}/contamination/pid_dedx')

    ROOT.gROOT.SetBatch(True)                     # Batch mode for ROOT
    ROOT.gStyle.SetErrorX(0)                      # Do not draw error along the X-axis
    ROOT.gStyle.SetOptStat(0)                     # Do not draw statistics box on the plots

    read_status('main')                           # Load file with the function trigger counters
    if os.path.exists(settings_file):             # Check and open the settings file
        config = open(settings_file)
    else:
        print(f'File not exist: {settings_file}') # Exit if the settings file does not exist
        exit()
    data = json.load(config)                      # Load data from the settings file

    Centrality = data['Events']['Centrality']     # Centrality bins
    n_centrality_bins = len(Centrality)           # Number of centrality bins
    Particles = data['Particles']                 # Particles information (PDG, mass, cuts -- everything)
    n_particles = len(Particles)                  # Number of particles
    PID_MODE = int(data['PID_mode'])              # PID mode
    print(f'PID mode: {PID_MODE}')

    INP_FILE = ROOT.TFile.Open(input_file,'READ') # Open and check (if not 'zombie') the input file
    check_root_file(INP_FILE)

    if(make_ptcorr):                              # Create file for the pT and pZ corrections
        corr_file = ROOT.TFile.Open(corrections_file,'RECREATE')
        check_root_file(corr_file)
        corr_file.cd()
    else:                                         # Open and check the efficiencies file
        EFF_FILE = ROOT.TFile.Open(efficiency_file,'READ')
        check_root_file(EFF_FILE)
        OUT_FILE = ROOT.TFile.Open(spectra_file, 'UPDATE') # Create and check the output file
        check_root_file(OUT_FILE)
        OUT_FILE.cd()

    c1 = TCanvas('c_landscape', '', 800, 600)     # 'Landscape' canvas
    c1.SetLeftMargin(0.14)                        # All four margins: left, right, bottom and top
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.14)
    c1.SetTopMargin(0.1)
    c2 = TCanvas('c_portrait', '', 600, 800)      # 'Portrait' canvas
    
    for particle in Particles:
        if(not make_ptcorr):
            os.system(f'mkdir -pv {outdir}/dca')
            os.system(f'mkdir -pv {outdir}/results/{{pt,dndy}}/{particle}') # Create the separate sub-directory
            os.system(f'mkdir -pv {outdir}/results/phasespace/{{uncorrected,corrected}}')
            os.system(f'mkdir -pv {outdir}/results/coal')                   # for the each particle in the file-system
            OUT_FILE.mkdir(f'{particle}', '', returnExistingDirectory=True) # and in the output file
            OUT_FILE.cd()
        for cbin in range(0, n_centrality_bins):
            common_suffix = f'{particle}_centrality{cbin}'
            if(make_ptcorr): # Copy pT and pZ TProfiles to the corrections file
                pz_corr_2d = INP_FILE.Get(f'pz__corr_{common_suffix}').Clone(f'p2d__pzcorr_{common_suffix}')
                pz_corr_2d.Write()
                pt_corr_2d = INP_FILE.Get(f'pt__corr_{common_suffix}').Clone(f'p2d__ptcorr_{common_suffix}')
                pt_corr_2d.Write()
            else:            # Make and analysis
                #draw_dca(particle, cbin, f'{outdir}/dca/{particle}_centrality{cbin}.{EXT}')
                #draw_dca_ratio(particle, cbin, f'{outdir}/dca/ratio_{particle}_centrality{cbin}.{EXT}')
                # --- TPC, Secondaries
                # ------- efficiency
                calculate_efficiency(f'h__eff_tpc_numerator_{particle}_centrality{cbin}',
                                     f'h__eff_tpc_denominator_{particle}_centrality{cbin}',
                                     f'h__efficiency_tpc_{particle}_centrality{cbin}', particle)
                # ------- contamination
                calculate_efficiency(f'h__cont_sec_numerator_{particle}_centrality{cbin}',
                                     f'h__cont_sec_denominator_{particle}_centrality{cbin}',
                                     f'h__contamination_secondaries_{particle}_centrality{cbin}', particle)
                # --- ToF
                # ------- efficiency
                calculate_efficiency(f'h__eff_pid_denominator_{particle}_centrality{cbin}',
                                     f'h__eff_tof_denominator_{particle}_centrality{cbin}',
                                     f'h__efficiency_tof_{particle}_centrality{cbin}', particle)
                if(PID_MODE == 0):
                    # --- evPID TPC only PID
                    # ------- efficiency
                    calculate_efficiency(f'h__eff_pid_tpc_numerator_{particle}_centrality{cbin}',
                                         f'h__eff_tof_denominator_{particle}_centrality{cbin}',
                                         f'h__efficiency_pid_tpc_{particle}_centrality{cbin}', particle)
                    # ------- contamination
                    calculate_efficiency(f'h__cont_pid_tpc_numerator_{particle}_centrality{cbin}',
                                         f'h__cont_pid_tpc_denominator_{particle}_centrality{cbin}',
                                         f'h__contamination_pid_tpc_{particle}_centrality{cbin}', particle)
                    # --- evPID ToF only PID
                    # ------- efficiency
                    calculate_efficiency(f'h__eff_pid_tof_numerator_{particle}_centrality{cbin}',
                                         f'h__eff_pid_denominator_{particle}_centrality{cbin}',
                                         f'h__efficiency_pid_tof_{particle}_centrality{cbin}', particle)
                    # ------- contamination
                    calculate_efficiency(f'h__cont_pid_tof_numerator_{particle}_centrality{cbin}',
                                         f'h__cont_pid_tof_denominator_{particle}_centrality{cbin}',
                                         f'h__contamination_pid_tof_{particle}_centrality{cbin}', particle)
                if(PID_MODE == 2):
                    # --- MpdPid combined PID
                    # ------- efficiency
                    eff = calculate_efficiency(f'h__cont_pid_denominator_{particle}_centrality{cbin}',
                                               f'h__eff_pid_denominator_{particle}_centrality{cbin}',
                                               f'h__efficiency_pid_{particle}_centrality{cbin}', particle)
                    # ------- Purity
                    eff = calculate_efficiency(f'h__eff_pid_numerator_{particle}_centrality{cbin}',
                                               f'h__cont_pid_denominator_{particle}_centrality{cbin}',
                                               f'h__purity_pid_{particle}_centrality{cbin}', particle)
                    # ------- contamination
                    eff = calculate_efficiency(f'h__cont_pid_numerator_{particle}_centrality{cbin}',
                                               f'h__cont_pid_denominator_{particle}_centrality{cbin}',
                                               f'h__contamination_pid_{particle}_centrality{cbin}', particle)
                    if DEDX:
                      # --- MpdPid TPC only (dE/dx) PID
                        # ------- efficiency
                        eff = calculate_efficiency(f'h__eff_pid_dedx_numerator_{particle}_centrality{cbin}',
                                                   f'h__eff_tof_denominator_{particle}_centrality{cbin}',
                                                   f'h__efficiency_pid_dedx_{particle}_centrality{cbin}', particle)
                        # ------- contamination
                        eff = calculate_efficiency(f'h__cont_pid_dedx_numerator_{particle}_centrality{cbin}',
                                                   f'h__cont_pid_dedx_denominator_{particle}_centrality{cbin}',
                                                   f'h__contamination_pid_dedx_{particle}_centrality{cbin}', particle)
                # MC spectra
                make_pt_spectra_mpdpid(0, 0, particle, cbin, f'h__pty_mc_{common_suffix}', 'mc')
                if(PID_MODE == 0):
                    # Corrected reconstructed spectra
                    make_pt_spectra_evpid('tpc', 0, particle, cbin, f'h__pty_pid_tpc_{common_suffix}') # Uncorrected
                    make_pt_spectra_evpid('tpc', 1, particle, cbin, f'h__pty_pid_tpc_{common_suffix}') # Corrected
                    make_pt_spectra_evpid('tof', 0, particle, cbin, f'h__pty_pid_tof_{common_suffix}') # Uncorrected
                    make_pt_spectra_evpid('tof', 1, particle, cbin, f'h__pty_pid_tof_{common_suffix}') # Corrected
                if(PID_MODE == 2):
                    # Uncorrected reconstructed spectra with combined PID = dE/dx and m^2
                    make_pt_spectra_mpdpid(0, 0, particle, cbin, f'h__pty_pid_{common_suffix}', 'uncorr')
                    # Corrected reconstructed spectra with combined PID = dE/dx and m^2
                    make_pt_spectra_mpdpid(1, 0, particle, cbin, f'h__pty_pid_{common_suffix}', 'corr')
                    # Combined efficiency: TPC + ToF + PID
                    make_overall_efficiency(0, particle, cbin)
                    if DEDX:
                        # Uncorrected reconstructed spectra with PID = dE/dx only
                        make_pt_spectra_mpdpid(0, 1, particle, cbin, f'h__pty_pid_dedx_{common_suffix}', 'uncorr')
                        # Corrected reconstructed spectra with PID = dE/dx only
                        make_pt_spectra_mpdpid(1, 1, particle, cbin, f'h__pty_pid_dedx_{common_suffix}', 'corr')                    
    INP_FILE.Close() # At this point the input file is not needed anymore
    if(make_ptcorr):
        corr_file.Close() # Same for the pT and pZ corrections case 
    else:
        EFF_FILE.Close()  # Efficiencies were already applied, efficiencies file can be closed
        for cbin in range(0,n_centrality_bins):
            cname = f'{Centrality[cbin][0]} - {Centrality[cbin][1]}%'
            for p_index, particle in enumerate(Particles):
                common_suffix = f'{particle}_centrality{cbin}'
                system = f'{SYS}, {particle}, {cname}'
                # Draw results: pT and rapidity spectra, efficiencies
                if(PID_MODE == 0):
                    draw_pt_spectra_evpid(c2, p_index, cbin, f'{outdir}/results/pt/{particle}/{particle}_centrality{cbin}')
            #         draw_phasespace(f'{particle}/phasespace_tpc_uncorr_{common_suffix}', system,
            #                         f'{outdir}/results/phasespace/uncorrected/h__pty_tpc_{common_suffix}.{EXT}')
            #         draw_phasespace(f'{particle}/phasespace_tpc_corr_{common_suffix}', system,
            #                         f'{outdir}/results/phasespace/corrected/h__pty_tpc_{common_suffix}.{EXT}')
            #         draw_phasespace(f'{particle}/phasespace_tof_uncorr_{common_suffix}', system,
            #                         f'{outdir}/results/phasespace/uncorrected/h__pty_tof_{common_suffix}.{EXT}')
            #         draw_phasespace(f'{particle}/phasespace_tof_corr_{common_suffix}', system,
            #                         f'{outdir}/results/phasespace/corrected/h__pty_tof_{common_suffix}.{EXT}')
            #     if(PID_MODE == 2):
            #         draw_phasespace(f'{particle}/phasespace_uncorr_{common_suffix}', system,
            #                         f'{outdir}/results/phasespace/uncorrected/h__pty_{common_suffix}.{EXT}')
            #         draw_phasespace(f'{particle}/phasespace_corr_{common_suffix}', system,
            #                         f'{outdir}/results/phasespace/corrected/h__pty_{common_suffix}.{EXT}')
            #         draw_pt_spectra_mpdpid(c2, 0, p_index, cbin, f'{outdir}/results/pt/{particle}/{particle}_centrality{cbin}')
            #         draw_overall_efficiency(c1, p_index, cbin, f'{outdir}/results/pt/{particle}/{common_suffix}_overall_efficiency')
            #         if particle != 'ap' and particle != 'He4':
            #             draw_dndy(0, c1, p_index, cbin, f'{outdir}/results/dndy/{particle}/{common_suffix}.{EXT}')
            #         if DEDX:
            #             draw_pt_spectra_mpdpid(c2, 1, p_index, cbin, f'{outdir}/results/pt/{particle}/{common_suffix}_dedx')
            #     if(extended):
            #         # TPC and primaries efficiency
            #         draw_efficiency(c1, f'{particle}/h__efficiency_tpc_{common_suffix}',
            #                             f'TPC efficiency, {system}',
            #                             f'{outdir}/efficiency/tpc/{common_suffix}.{EXT}')
            #         draw_efficiency_midrapidity(c1, f'{particle}/h__efficiency_tpc_{common_suffix}',
            #                                         f'TPC efficiency at |y| < 0.5, {system}',
            #                                         f'{outdir}/efficiency/tpc/mid_{common_suffix}.{EXT}')
            #         # Secondaries contamination
            #         draw_efficiency(c1, f'{particle}/h__contamination_secondaries_{common_suffix}',
            #                             f'Secondaries contamination, {system}',
            #                             f'{outdir}/contamination/secondaries/{common_suffix}.{EXT}')
            #         draw_efficiency_midrapidity(c1, f'{particle}/h__contamination_secondaries_{common_suffix}',
            #                                         f'Secondaries contamination at |y| < 0.5, {system}',
            #                                         f'{outdir}/contamination/secondaries/mid_{common_suffix}.{EXT}')
            #         # ToF efficiency
            #         draw_efficiency(c1, f'{particle}/h__efficiency_tof_{common_suffix}',
            #                             f'ToF efficiency, {system}',
            #                             f'{outdir}/efficiency/tof/{common_suffix}.{EXT}')
            #         draw_efficiency_midrapidity(c1, f'{particle}/h__efficiency_tof_{common_suffix}',
            #                                         f'ToF efficiency at |y| < 0.5, {system}',
            #                                         f'{outdir}/efficiency/tof/mid_{common_suffix}.{EXT}')
            #         if(PID_MODE == 0):
            #             # ecPID TPC only efficiency
            #             draw_efficiency(c1, f'{particle}/h__efficiency_pid_tpc_{common_suffix}',
            #                                 f'evPID efficiency (TPC only), {system}',
            #                                 f'{outdir}/efficiency/pid/tpc_{common_suffix}.{EXT}')
            #             # ecPID TPC only contamination
            #             draw_efficiency(c1, f'{particle}/h__contamination_pid_tpc_{common_suffix}',
            #                                 f'evPID contamination (TPC only, {system}',
            #                                 f'{outdir}/contamination/pid/tpc_{common_suffix}.{EXT}')
            #             # ecPID ToF only efficiency
            #             draw_efficiency(c1, f'{particle}/h__efficiency_pid_tof_{common_suffix}',
            #                                 f'evPID efficiency (ToF only), {system}',
            #                                 f'{outdir}/efficiency/pid/tof_{common_suffix}.{EXT}')
            #             # ecPID ToF only contamination
            #             draw_efficiency(c1, f'{particle}/h__contamination_pid_tof_{common_suffix}',
            #                                 f'evPID contamination (ToF only, {system}',
            #                                 f'{outdir}/contamination/pid/tof_{common_suffix}.{EXT}')
            #         if(PID_MODE == 2):
            #             # MpdPid combined efficiency
            #             draw_efficiency(c1, f'{particle}/h__efficiency_pid_{common_suffix}',
            #                                 f'PID efficiency, {system}',
            #                                 f'{outdir}/efficiency/pid/{common_suffix}.{EXT}', 1)
            #             draw_efficiency_midrapidity(c1, f'{particle}/h__efficiency_pid_{common_suffix}',
            #                                             f'PID efficiency at |y| < 0.5, {system}',
            #                                             f'{outdir}/efficiency/pid/mid_{common_suffix}.{EXT}')
            #             # MpdPid combined purity
            #             draw_efficiency(c1, f'{particle}/h__purity_pid_{common_suffix}',
            #                                 f'PID purity, {system}',
            #                                 f'{outdir}/efficiency/pid/purity_{common_suffix}.{EXT}')
            #             draw_efficiency_midrapidity(c1, f'{particle}/h__purity_pid_{common_suffix}',
            #                                             f'PID purity at |y| < 0.5, {system}',
            #                                             f'{outdir}/efficiency/pid/mid_purity_{common_suffix}.{EXT}')
            #             # MpdPid combined contamination
            #             draw_efficiency(c1, f'{particle}/h__contamination_pid_{common_suffix}',
            #                                 f'PID contamination, {system}',
            #                                 f'{outdir}/contamination/pid/{common_suffix}.{EXT}')
            #             draw_efficiency_midrapidity(c1, f'{particle}/h__contamination_pid_{common_suffix}',
            #                                             f'PID contamination at |y| < 0.5, {system}',
            #                                             f'{outdir}/contamination/pid/mid_{common_suffix}.{EXT}')
            #             if DEDX:
            #                 # MpdPid dE/dx efficiency
            #                 draw_efficiency(c1, f'{particle}/h__efficiency_pid_dedx_{common_suffix}',
            #                                     f'PID (dE/dx only) efficiency, {system}',
            #                                     f'{outdir}/efficiency/pid_dedx/{common_suffix}.{EXT}')
            #                 # MpdPid dE/dx contamination
            #                 draw_efficiency(c1, f'{particle}/h__contamination_pid_dedx_{common_suffix}',
            #                                     f'PID (dE/dx only) contamination, {system}',
            #                                     f'{outdir}/contamination/pid_dedx/{common_suffix}.{EXT}')
            # if(PID_MODE == 2):
            #     draw_b2(c1, cbin, f'{outdir}/results/coal/b2_centrality{cbin}')        # Calculate and draw the B2 coalescence parameter
            #    	draw_b3he3(c1, cbin, f'{outdir}/results/coal/b3he3_centrality{cbin}')  # Calculate and draw the B3 (for 3He) coalescence parameter
        OUT_FILE.Close()
    config.close()
    if not make_ptcorr and report_file is not None: # Produce the report file (both tex and pdf)
        generate_report(extended, outdir, Path(report_file).stem)



##  This subroutine draws the efficiency or contamination histograms
#  \param canvas TCanvas prepared for drawing
#  \param psname The name of the selected efficiency (or contamination) histogram in the input file
#  \param title The title for the output histogram
#  \param fname The output pdf-file name
#  \param maximum [optional] The histogram maximum for the Z-axis
def draw_efficiency(canvas, psname, title, fname, maximum = None):
    order = 3                                                 # These three lines are common for several functions
    CURRENT_ST[order] += 1                                    # Increment the function calls counts
    function_name = sys._getframe().f_code.co_name            # Here we get the function name, the execution order and a counter
    if skip_function(function_name, order): return            # Thus we can skip this function execution if it was already
                                                              #   called during the previous program run
    canvas.cd()                                               # Go to this canvas
    hPhaseSpace = OUT_FILE.Get(psname).Clone('hPhaseSpace')   # Get the needed phase space to draw
    frame = canvas.DrawFrame(-2.99, 0, 2.99, 4.99)            # Draw frame with chosen X/Y axes limits
    make_fancy_frame(frame, title, 'y', 'p_{T}, GeV/c')       # Set frame title, axes labels and fonts
    if maximum:                                               # If the maximum parameter is set
        hPhaseSpace.SetMaximum(maximum)                       #   then re-define the histogram maximum of the Z-axis
    hPhaseSpace.Draw('same colz')                             # Draw phase-space inside the frame
    canvas.Print(fname)                                       # The canvas content is printed to the file
    write_status(f'{function_name} {CURRENT_ST[order]}')      # Here the function execution counter is updated in the status file



##  This subroutine draws the efficiency or contamination histograms at midrapidity
#  \param canvas TCanvas prepared for drawing
#  \param psname The name of the selected efficiency (or contamination) histogram in the input file
#  \param title The title for the output histogram
#  \param fname The output pdf-file name
def draw_efficiency_midrapidity(canvas, psname, title, fname):
    order = 11
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    rb_low  = -0.5                            # It's a midrapidity function
    rb_high =  0.5                            #    so, rapidity limits are hard-coded
    rb_width = math.fabs(rb_high - rb_low)    # But lets calculate the midrapidity bin width :)
    canvas.cd()
    hPhaseSpace = OUT_FILE.Get(psname).Clone('hPhaseSpace')                               # Here we select the needed phase-space
    for hbin in range(1, hPhaseSpace.GetNbinsX()):                                        # In this loop we find 
        if(math.fabs(hPhaseSpace.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA): #   the number of the left rapidity bin:
            left_edge = hbin                                                              #   so  'left_edge'  is the bin number!
            break
    right_edge = int(left_edge + rb_width / hPhaseSpace.GetXaxis().GetBinWidth(left_edge) - 1) # Then we calculate the bin number for the right limit
    hMidEff = hPhaseSpace.ProjectionY('midrapidity_efficiency', left_edge, right_edge)         # This is how the projection of the phase-space on the Y-axis is made
    hMidEff.Scale(1. / (right_edge - left_edge + 1))                                           # And now we scale the resulting projection by the number of the
                                                                                               #   projected bins (otherwise it will be the sum over all these bins).
    frame = canvas.DrawFrame(0, 0, 3.99, 1.09)                      # Again some fancy stuff
    make_fancy_frame(frame, title, 'p_{T}, GeV/c', 'Fraction')      #   for the frame
    make_fancy_histogram(hMidEff, MS_MC, LC_MC, 0, 0, msize = 2)    #   and for the histogram
    hMidEff.Draw('same hist P')
    canvas.Print(fname)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the phase-space histograms (very similar to the draw_efficiency() function)
#  \param psname The name of the selected phase-space
#  \param title The title for the output histogram
#  \param fname The output pdf-file name
#  \param minimum [optional] The histogram minimum for the Z-axis
#  \param maximum [optional] The histogram maximum for the Z-axis
def draw_phasespace(psname, title, fname, minimum = None, maximum = None):
    order = 7                                                 # These three lines are common for several functions
    CURRENT_ST[order] += 1                                    # Increment the function calls counts
    function_name = sys._getframe().f_code.co_name            # Here we get the function name, the execution order and a counter
    if skip_function(function_name, order): return            # Thus we can skip this function execution if it was already
                                                              # called during the previous program run
    canvas = TCanvas()                                        # The new TCanvas is defined
    canvas.SetLeftMargin(0.14)                                #   with some custom margins
    canvas.SetRightMargin(0.1)
    canvas.SetBottomMargin(0.14)
    canvas.SetTopMargin(0.1)
    canvas.SetLogz(1)
    canvas.cd()                                               # Go to this canvas
    hPhaseSpace = OUT_FILE.Get(psname).Clone('hPhaseSpace')   # Get the needed phase space to draw
    frame = canvas.DrawFrame(-2.99, 0, 2.99, 4.99)            # Draw frame with chosen X/Y axes limits
    make_fancy_frame(frame, title, 'y', 'p_{T}, GeV/c')       # Set frame title, axes labels and fonts
    if minimum:                                               # If 'minimum' is defined
        hPhaseSpace.SetMinimum(minimum)                       #   use for the histogram Z-axis
    if maximum:                                               # If 'maximum' is defined
        hPhaseSpace.SetMaximum(maximum)                       #   use for the histogram Z-axis maximum
    hPhaseSpace.Draw('same colz')                             # Draw phase-space inside the frame
    canvas.Print(fname)                                       # The canvas content is printed to the file
    write_status(f'{function_name} {CURRENT_ST[order]}')      # Here the function execution counter is updated in the status file



##  This subroutine draws the DCA distribution
#  \param particle Selected particle
#  \param cbin Selected centrality bin
#  \param fname The output pdf-file name
def draw_dca(particle, cbin, fname):
    order = 12
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    canvas = TCanvas()                                        # The new TCanvas is defined
    canvas.SetLeftMargin(0.14)                                #   with some custom margins
    canvas.SetRightMargin(0.1)
    canvas.SetBottomMargin(0.14)
    canvas.SetTopMargin(0.1)
    canvas.SetLogy(1)
    canvas.cd()
    hDCAsec  = INP_FILE.Get(f'h__dca_sec_{particle}_centrality{cbin}').Clone('hDCAsec')   # DCA of the secondary particles (by the MC info)
    hDCAprim = INP_FILE.Get(f'h__dca_prim_{particle}_centrality{cbin}').Clone('hDCAprim') # DCA of the primary particles (by the MC info)
    hDCA     = INP_FILE.Get(f'h__dca_{particle}_centrality{cbin}').Clone('hDCA')          # DCA for the accepted by the track quality cuts ('primary') particles
    hDCAsum  = hDCAprim.Clone('hDCAsum')
    hDCAsum.Add(hDCAsec)                                                                  # Sum of the secondaries and primaries
    frame = canvas.DrawFrame(0, 1E1, 6.99, hDCAsum.GetMaximum()*1E1)
    make_fancy_frame(frame, f'DCA for {particle}', 'DCA, cm', 'counts')
    make_fancy_histogram(hDCAsum,  0, 0, 1, 1, lwidth = 3) # Markers, lines and colours for histograms
    make_fancy_histogram(hDCAprim, 0, 0, 2, 3, lwidth = 3) # Markers, lines and colours for histograms
    make_fancy_histogram(hDCAsec,  0, 0, 2, 4, lwidth = 3) # Markers, lines and colours for histograms
    make_fancy_histogram(hDCA,     0, 0, 3, 2, lwidth = 3) # Markers, lines and colours for histograms
    hDCAsum.Draw('same hist')
    hDCAprim.Draw('same hist')
    hDCAsec.Draw('same hist')
    hDCA.Draw('same hist')
    ratio_sel    = hDCA.Integral() / hDCAsum.Integral() * 100
    legend = TLegend(0.5, 0.7, 0.88, 0.9)    # Create, modify and place legend
    legend.SetNColumns(1)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.AddEntry(hDCAsum,  'Total', 'L')
    legend.AddEntry(hDCAprim, 'Primary', 'L')
    legend.AddEntry(hDCAsec,  'Secondary', 'L')
    legend.AddEntry(hDCA,     'Accepted by DCA ({:.1f}%)'.format(ratio_sel), 'L')
    legend.Draw()
    tbox = ROOT.TPaveText(0.15, 0.05, 0.45, 0.5, 'NDC')       # The text box can be added to the plot. It contains:
    tbox.AddText('Total = {:.3f}'.format(hDCAsum.Integral())) # The sum of the primaries and secondaries
    tbox.AddText('Prim = {:.3f}'.format(hDCAprim.Integral())) # The number of the primary particles
    tbox.AddText('Sec = {:.3f}'.format(hDCAsec.Integral()))   # The number of the secondary particles
    tbox.AddText('DCA = {:.3f}'.format(hDCA.Integral()))      # The number of the particles accepted by the DCA track quality cut
    # ~ tbox.Draw()
    canvas.Print(fname)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the ratio of the DCA distributions
#  \param particle Selected particle
#  \param cbin Selected centrality bin
#  \param fname The output pdf-file name
def draw_dca_ratio(particle, cbin, fname):
    order = 13
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    canvas = TCanvas()
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.1)
    canvas.SetBottomMargin(0.14)
    canvas.SetTopMargin(0.1)
    canvas.SetLogy(0)
    canvas.cd()
    hDCAsec  = INP_FILE.Get(f'h__dca_sec_{particle}_centrality{cbin}').Clone('hDCAsec')   # DCA of the secondary particles (by the MC info)
    hDCAprim = INP_FILE.Get(f'h__dca_prim_{particle}_centrality{cbin}').Clone('hDCAprim') # DCA of the primary particles (by the MC info)
    hDCA     = INP_FILE.Get(f'h__dca_{particle}_centrality{cbin}').Clone('hDCA')          # DCA for the accepted by the track quality cuts ('primary') particles
    hDCAsum  = hDCAprim.Clone('hDCAsum')
    hDCAsum.Add(hDCAsec)                                                                  # Sum of the secondaries and primaries

    hDCAsec.Divide(hDCAsum)  # Ratio: secondaries / (secondaries + primaries)
    hDCAprim.Divide(hDCAsum) # Ratio: primaries /  (secondaries + primaries)

    frame = canvas.DrawFrame(0, 0, 6.99, 1.09)
    make_fancy_frame(frame, f'DCA for {particle}', 'DCA, cm', 'ratio')
    make_fancy_histogram(hDCAsum,  0, 0, 1, 1, lwidth = 3) # Markers, lines and colours for histograms
    make_fancy_histogram(hDCAprim, 0, 0, 2, 3, lwidth = 3) # Markers, lines and colours for histograms
    make_fancy_histogram(hDCAsec,  0, 0, 2, 4, lwidth = 3) # Markers, lines and colours for histograms
    make_fancy_histogram(hDCA,     0, 0, 1, 2, lwidth = 3) # Markers, lines and colours for histograms
    hDCAprim.Draw('same hist')
    hDCAsec.Draw('same hist')
    legend = TLegend(0.5, 0.7, 0.9, 0.9)    # Create, modify and place legend
    legend.SetNColumns(1)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.AddEntry(hDCAprim, 'Primary / Total', 'L')
    legend.AddEntry(hDCAsec,  'Secondary / Total', 'L')
    legend.Draw()
    canvas.Print(fname)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine is a thermal fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum \f$ p_{T} \f$
#  \param par The fit parameters set:
#     - par[0] dN/dy
#     - par[1] Temperature T
#     - par[2] Particle mass (must be fixed)
#  \return The fit value
def fit_thermal(x, par):
    dndy = par[0]
    T    = par[1]
    m0   = par[2]
    mt   = math.sqrt(x[0]*x[0] + m0*m0)
    val  = (dndy/(T*(m0 + T))) * math.exp(-(mt - m0)/T)
    return val



##  This subroutine is a Blast-Wave fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum \f$ p_{T} \f$
#  \param par The fit parameters set:
#     - par[0] Blastwave normalization constant
#     - par[1] Mass of the particle (must be fixed)
#     - par[2] Blastwave temperature
#     - par[3] Blastwave flow velocity
#  \return The fit value
def fit_blastwave_pt_invariant(x, par):
    pt    = x[0]
    C     = par[0]
    m0    = par[1]
    T     = par[2]
    beta  = par[3]
    mt    = ROOT.TMath.Sqrt(pt*pt + m0*m0)
    R_max = 20
    steps = 40
    dr    = R_max / steps

    fitval = 0
    for i in range(steps):
        r      = i * dr
        betar  = beta * ROOT.TMath.Power(r/R_max, 1.0)
        rho    = TMath.ATanH(betar)
        b1     = ROOT.TMath.BesselI0( pt * ROOT.TMath.SinH(rho) / T )
        b2     = ROOT.TMath.BesselK1( mt * ROOT.TMath.CosH(rho) / T )
        fitval += r * mt * dr * b1 * b2
    return C * fitval



##  This subroutine is a Blast-Wave fit for the \f$ d^{2}N/dp_{T}dy \f$ spectra
#  \param x Transverse momentum $p_{T}$
#  \param par The fit parameters set:
#     - par[0] Blastwave normalization constant
#     - par[1] Mass of the particle (must be fixed)
#     - par[2] Blastwave temperature
#     - par[3] Blastwave flow velocity
#  \return The fit value
def fit_blastwave_pt(x, par):
    pt    = x[0]
    C     = par[0]
    m0    = par[1]
    T     = par[2]
    beta  = par[3]

    mt    = ROOT.TMath.Sqrt(pt*pt + m0*m0)
    R_max = 20
    steps = 40
    dr    = R_max / steps

    fitval = 0
    for i in range(steps):
        r      = i * dr
        betar  = beta * ROOT.TMath.Power(r/R_max, 1.0)
        rho    = TMath.ATanH(betar)
        b1     = ROOT.TMath.BesselI0( pt * ROOT.TMath.SinH(rho) / T )
        b2     = ROOT.TMath.BesselK1( mt * ROOT.TMath.CosH(rho) / T )
        fitval += r * mt * dr * b1 * b2
    return C * pt * fitval



##  This subroutine is a double-exponential (thermal) fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum $p_{T}$
#  \param par The fit parameters set:
#     - par[0] Constant 1
#     - par[1] Constant 2
#     - par[2] Temperature T 1
#     - par[3] Temperature T 2
#     - par[4] Particle mass (must be fixed)
#  \return The fit value
def fit_2exp(x, par):
    pt  = x[0]   # Transverse momentum
    m0  = par[4] # Particle mass (fixed parameter)
    mtm = ROOT.TMath.Sqrt(pt * pt + m0 * m0) - m0
    try:
        x1 = (1. + par[1])/(par[2] * (m0 + par[2])) / ROOT.TMath.Exp(mtm / par[2]) # This can be zero
    except ZeroDivisionError:                                                      # So, the exception is made
        x1 = 0
    try:
        x2 = par[1] / (par[3] * (m0 + par[3])) / ROOT.TMath.Exp(mtm / par[3])      # This too can be zero
    except ZeroDivisionError:                                                      # So, the exception is made too
        x2 = 0
    val = par[0] * (x1 - x2)
    return val



##  This subroutine is a Bose-Einstein fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum \f$ p_{T} \f$
#  \param par The fit parameters set:
#     - par[0] dN/dy
#     - par[1] Temperature T
#     - par[2] Particle mass (must be fixed)
#  \return The fit value
def fit_bose_einstein(x, par):
    pt   = x[0]
    dndy = par[0]
    T    = par[1]
    m0   = par[2]
    mt   = math.sqrt(x[0]*x[0] + m0*m0)
    val  = dndy / (math.exp(mt/T) - 1)
    return val



## The  wrapper class for the TF1 fit functions
#
# It allows to apply fits simultaneously using the multiprocessing.Pool.
# See the explanations here:
# https://root-forum.cern.ch/t/multiprocessing-fits-within-pyroot/60297/10  
# https://github.com/root-project/root/issues/16184
class TF1Wrapper:
    def __init__(self, *args):
        self.f = ROOT.TF1(*args)
        self.args = args



##  This subroutine applies the fit to the histogram
#  \param hist TH1 histogram
#  \param fit TF1 fit function
#  \param frange The fit range ([low, high])
#  \return The TF1 fit
def apply_fit(hist, fit, *frange):
    fitf = ROOT.TF1(*fit.args)          # "Workaround": we are unfortunately forced to recreate the function
    for npar in range(fit.f.GetNpar()):
        fitf.SetParameter(npar, fit.f.GetParameter(npar))
    fitname = fitf.GetName()            # The fit name
    # Fix the parameters accordingly to the fit name
    if fitname == 'blastwave_inv':
        fitf.FixParameter(1, fit.f.GetParameter(1)) # Particle mass must be always fixed
        fitf.SetParLimits(2, 0.001, 0.99)           # Blastwave temperature
        fitf.SetParLimits(3, 0.001, 0.99)           # Blastwave flow velocity
    elif fitname == 'thermal':
        fitf.FixParameter(2, fit.f.GetParameter(2)) # Particle mass must be always fixed
        fitf.SetParLimits(1, 0.001, 1E3)            # The temperature parameter T should not be negative
    elif fitname == 'fit_dexp':
        fitf.FixParameter(4, fit.f.GetParameter(4)) # Particle mass must be always fixed
    elif fitname == 'bosein':
        fitf.FixParameter(2, fit.f.GetParameter(2)) # Particle mass must be always fixed
        fitf.SetParLimits(1, 0.001, 1E3)            # The temperature parameter T should not be negative

    counter = 0 # The fit tries counter
    while True:
        try:
            res = hist.Fit(fitf, '0Q', '', *frange)                     # Let's make the first attempt of the fit
        except ZeroDivisionError:                                       # If the 'ZeroDivisionError' occurs
            res = -1                                                    #   Set the fit status to the '-1'
            f0 = fitf.GetParameter(0)                                   #   Get the first (actually, the '0') parameter 
            fitf.SetParameter(0, f0 + random.uniform(-f0*0.9, f0*0.9))  #     and randomly modify it slightly within the +/- range of its own value and perform the fit again
        if (int(res) == 3):                                             # If the fit result is equal to '3', that mean some error occurs and the fit is not good, so,
            f0 = fitf.GetParameter(0)                                   #   again, get the first (actually, the '0') parameter 
            fitf.SetParameter(0, f0 + random.uniform(-f0*0.9, f0*0.9))  #     and randomly modify it slightly within the +/- range of its own value and perform the fit again
        if (int(res) >= 0 and int(res) != 3):                           # Now, if the fit result is '0' or some positive integer except '3'
            break                                                       #   the fit procedure is done -- exit from the fit loop
        counter += 1                                                    # Increase the fit attempts counter
        if counter > 15:                                                # If there was 15 attemtps to fix or to improve the fit
            break                                                       #   exit from the fit loop
    print(f'name={fitf.GetName()}, status={int(res)}, try={int(counter)}')
    return fitf



##  This subroutine draws the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra for the defined rapidity bins
#  \param canvas TCanvas prepared for drawing
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
def draw_pt_spectra_mpdpid(canvas, do_dedx, p_pos, cbin, fname):
    order = 1
    function_name = sys._getframe().f_code.co_name
    CURRENT_ST[order] += 1
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]                         # Particle name
    cbins    = list(Centrality)[cbin]                         # Centrality bin limits
    cname    = f'{cbins[0]} - {cbins[1]}%'                    # Centrality bin name
    pt_start = list(Particles.values())[p_pos]['pt_bins'][1]  # pT limits: low and 
    pt_end   = list(Particles.values())[p_pos]['pt_bins'][2]  #   high
    p_mass   = float(list(Particles.values())[p_pos]['Mass']) # Particle mass
    print(f' draw pT =============== {pname},  {cname},  {pt_start} < pT < {pt_end}  =============== ')
    
    if(do_dedx):
        postfix='_dedx'                                       # Histograms name postfix for the 'PID = dE/dx only'
    else:
        postfix=''
    canvas.cd()
    vMC = []                                                  # The list of Monte-Carlo histograms
    vUncorr = []                                              # The list of uncorrected histograms
    vCorr = []                                                # The list of corrected histograms
    vDummy = []                                               # Stupid hack for drawing
    common_suffix = f'{pname}_centrality{cbin}'               # Histograms name suffix for each particle and centrality
    for rbin in range(len(rapidity_bins)):                    # Populate lists with histograms for the eah rapidity bin
        common_postfix = f'{postfix}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hMC     = OUT_FILE.Get(f'{pname}/h__pt_{common_suffix}_mc_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}').Clone('hMC')
        hUncorr = OUT_FILE.Get(f'{pname}/h__pt_{common_suffix}_uncorr{common_postfix}').Clone('hUncorr')
        hCorr   = OUT_FILE.Get(f'{pname}/h__pt_{common_suffix}_corr{common_postfix}').Clone('hCorr')
        hDummy  = hCorr.Clone('hDummy')
        vMC.append(hMC)
        vUncorr.append(hUncorr)
        vCorr.append(hCorr)
        vDummy.append(hDummy)

    # Blast-Wave (BW) fit for the pT invariant spectra
    blastwave_inv = TF1Wrapper('blastwave_inv', fit_blastwave_pt_invariant, pt_start, pt_end, 4)
    blastwave_inv.f.SetParNames('C', 'm0', 'T', 'beta')
    blastwave_inv.f.SetParameters(100, 0.938, 0.15, 0.6)
    blastwave_inv.f.SetParLimits(2, 0.001, 0.99)
    blastwave_inv.f.SetParLimits(3, 0.001, 0.99)

    # Blast-Wave fit for the pT spectra
    blastwave = ROOT.TF1('blastwave', fit_blastwave_pt, pt_start, pt_end, 4)

    # Thermal fit for the pT invariant spectra
    thermal = TF1Wrapper('thermal', fit_thermal, pt_start, pt_end, 3)
    thermal.f.SetParNames('dN/dy', 'T', 'm0')
    thermal.f.SetParameters(20, 0.1, 0.938)

    # Double exponential fit for the pT invariant spectra
    fit_dexp = TF1Wrapper('fit_dexp', fit_2exp, pt_start, pt_end, 5)
    fit_dexp.f.SetParameters(30., 0.4, 0.2, 0.1, 0.938)


    # Bose-Einstein fit for the pT invariant spectra
    bosein = TF1Wrapper('bosein', fit_bose_einstein, pt_start, pt_end, 3)
    bosein.f.SetParNames('dN/dy', 'T', 'm0')
    bosein.f.SetParameters(40, 0.1, 0.938)

    # Fit ranges and some initial conditions
    plot_range = [0, 6]                                       # Common histograms plotting range
    if (pname == 'pim' or pname == 'pip'):   # pi- and pi+ fits parameters set
        plot_range      = [0, 1.59]                           # Plot range for the selected particle
        fit_range       = [0.1, 2.0]                          # BW fit range for the selected particle
        fit_range_dexp  = fit_range                           # Double exponential fit range for the selected particle
        dndy_exp_limits = [0.1, 1.5]                          # Range for the dN/dy calculation using the corrected pT spectra
        blastwave_inv.f.SetParameters(100, p_mass, 0.2, 0.2)  # BW fit parameters for the selected particle
        thermal.f.SetParameters(20, 0.1, p_mass)              # Thermal fit parameters for the selected particle
        fit_dexp.f.SetParameters(50, 5, 0.3, 0.3, p_mass)     # Double exponential fit parameters for the selected particle
    elif (pname == 'km' or pname == 'kp'):   # K- and K+ -- same as above
        plot_range      = [0, 1.69]                    
        fit_range       = [0.2, 1.5]
        fit_range_dexp  = fit_range
        dndy_exp_limits = [0.2, 1.5]
        blastwave_inv.f.SetParameters(100, p_mass, 0.2, 0.2)
        fit_dexp.f.SetParameters(20, 2, 0.3, 0.3, p_mass)
    elif (pname == 'p'):                     # p -- same as above
        plot_range      = [0, 2.19]
        fit_range       = [0.5, 1.5]
        fit_range_dexp  = fit_range
        dndy_exp_limits = [0.5, 1.5]
        blastwave_inv.f.SetParameters(1000, p_mass, 0.1, 0.1)
    elif (pname == 'd'):                     # d -- same as above
        plot_range      = [0, 3.09]
        fit_range       = [0.9, 2.5]
        fit_range_dexp  = [1.0, 4.0]
        dndy_exp_limits = [0.8, 2.6]
        blastwave_inv.f.SetParameters(1000, p_mass, 0.2, 0.2)
        fit_dexp.f.SetParameters(1, 24, 0.26, 0.22, p_mass)
    elif (pname == 't'):                     # t -- same as above
        plot_range      = [0, 3.99]
        fit_range       = [1.4, 3.5]
        fit_range_dexp  = [1.5, 4.0]
        dndy_exp_limits = [1.2, 3.9]
        thermal.f.SetParameters(10, 0.1, p_mass)
        bosein.f.SetParameters(10, 0.2, p_mass)
        fit_dexp.f.SetParameters(0.2, 20, 0.4, 0.4, p_mass)
    elif (pname == 'He3'):                    # 3He -- same as above
        plot_range      = [0, 3.99]
        fit_range       = [1.4, 3.5]
        fit_range_dexp  = [1.5, 4.0]
        dndy_exp_limits = [1.2, 3.9]
        thermal.f.SetParameters(5, 0.1, p_mass)
        bosein.f.SetParameters(5, 0.2, p_mass)
        fit_dexp.f.SetParameters(0.2, 20, 0.4, 0.4, p_mass)

    # Set the right particle mass for the each fit
    blastwave_inv.f.FixParameter(1, p_mass)
    thermal.f.FixParameter(2, p_mass)
    fit_dexp.f.FixParameter(4, p_mass)
    bosein.f.FixParameter(2, p_mass)

    canvas.SetLogy(1)                       # Set the canvas Y-axis as log
    legend = TLegend(0.7, 0.6, 0.9, 0.9)    # Create, modify and place legend
    legend.SetNColumns(1)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)

    for rbin in range(len(rapidity_bins)):
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'  # Output file name for the each plot
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'

        ratio_cor = TRatioPlot(vDummy[rbin], vMC[rbin])                              # Corrected ratio plot
        ratio_uncor = TRatioPlot(vUncorr[rbin], vMC[rbin])                           # Uncorrected ratio plot
        ratio_cor.SetSeparationMargin(0.0)                                           # The margin between upper and lower (ratio) part
        ratio_cor.SetRightMargin(0.01)                                               # The right margin
        ratio_cor.SetLeftMargin(0.12)                                                # The left margin

        make_fancy_histogram(vMC[rbin], MS_MC, LC_MC, 1, LC_MC)                      # Markers, lines and colours for histograms:
        make_fancy_histogram(vUncorr[rbin], MS_UNCOR, LC_UNCOR, 1, LC_UNCOR)         # MC, uncorrected and corrected
        make_fancy_histogram(vCorr[rbin], MS_COR, LC_COR, 1, LC_COR)
        make_fancy_histogram(vDummy[rbin], 0, 0, 0, 0)

        if(rbin == 0):                                                               # The legend must be created one
            legend.AddEntry(vMC[rbin], 'MC', 'P')                                    #   but then it can be used for the each rapidity bin
            # ~ legend.AddEntry(vUncorr[rbin], 'Uncorrected', 'PL')
            legend.AddEntry(vCorr[rbin], 'Corrected', 'PL')
            if((pname != 'ap' and pname != 'He4') and postfix==''):                  # Here the anti-protons and the He4
                obj = ROOT.TGraph()                                                  #   are not analyzed
                make_fancy_histogram(obj, 23,  1, 2, 14, 0.5, 3)
                legend.AddEntry(obj, 'Blast-Wave', 'PL')                             # Add BW fit to the legend
                objt = ROOT.TGraph()
                make_fancy_histogram(objt, 25,  6, 3, 6, 0.5, 3)
                legend.AddEntry(objt, 'Thermal', 'PL')                               # Add thermal fit to the legend
                obje = ROOT.TGraph()
                make_fancy_histogram(obje, 24,  4, 3, 4, 0.5, 3)
                legend.AddEntry(obje, 'Double exp.', 'PL')                           # Add double-exponential fit to the legend
                objb = ROOT.TGraph()
                make_fancy_histogram(objb, 26,  42, 3,  42, 0.5, 3)
                legend.AddEntry(objb, 'Bose-Einstein', 'PL')                         # Add double-exponential fit to the legend

        ratio_uncor.Draw()                                      # Draw uncorrected plot to have this object in the future
        g_uncor = ratio_uncor.GetLowerRefGraph()                # Get the lower (ratio) part of the uncorrected reco results
        g_uncor.SetLineColor(LC_UNCOR)                          # Set the uncorrected results line colour
        ratio_cor.SetH1DrawOpt('0')                             
        ratio_cor.Draw('')                                      # Draw corrected results, all futher histograms will be
        ROOT.gPad.Update()                                      #   plotted on these pads (lower and upper)
        ratio_cor.GetLowerRefGraph().SetLineColor(LC_COR)       # Set the corrected results line colour
        ratio_cor.GetLowerRefGraph().SetLineWidth(2)            # and the line width
        ROOT.gPad.Update()
        y_low = ratio_cor.GetUpperPad().GetFrame().GetY1()      # The upper pad Y-axis maximum
        y_max = ratio_cor.GetUpperPad().GetFrame().GetY2()      #   and minimum
        ratio_cor.GetLowerRefGraph().SetMinimum(0)              # Set the lower pad minimum to 0
        ratio_cor.GetLowerRefGraph().SetMaximum(2)              # Set the lower par maximum to 2
        ratio_cor.GetLowerRefXaxis().SetRangeUser(*plot_range)  # Set the lower (and tus upper too) pad the defined above plot range
        if(p_pos == 7):                                  # Special case for tritons
            ratio_cor.GetLowerRefGraph().SetMaximum(4)
        if(p_pos >= 6):                                  # Special case for d and above
            ratio_cor.GetUpperRefYaxis().SetRangeUser(math.pow(10, y_max - 3), math.pow(10, y_max+0.2)) # 3 orders of magnitude less
        if(p_pos >= 8):                                  # Special case for He3 and He4
            ratio_cor.GetUpperRefYaxis().SetRangeUser(math.pow(10, y_max - 3.5), math.pow(10, y_max+0.5)) # 3.5 orders of magnitude less
        ratio_cor.GetLowerRefYaxis().SetNdivisions(6)
        ratio_cor.GetUpperPad().cd()                            # Go to the upper pad and draw histograms and fits there
        vMC[rbin].Draw('same hist PE')                          # Monte-Carlo histogram for the current rapidity bin
        # ~ vUncorr[rbin].Draw('same hist PE')                      # Uncorrected results histogram for the current rapidity bin
        for hbin in range(1, vCorr[rbin].GetSize() - 1):
            content_mc = vMC[rbin].GetBinContent(hbin)
            if content_mc == 0:
                continue
            content_rc = vCorr[rbin].GetBinContent(hbin)
            ratio_rcmc = content_rc / content_mc
            if (ratio_rcmc < 0.8 or ratio_rcmc > 1.2):          # If the difference between the reconstructed and MC points is high:
                vCorr[rbin].SetBinContent(hbin, 0)              #   do not draw these reconstructed points
        vCorr[rbin].Draw('same hist PE')                        # Corrected results histogram for the current rapidity bin
        vCorr[rbin].Draw('same hist')                           # Another stupid hack for drawing
        if((pname != 'ap' and pname != 'He4') and postfix==''):
            pool = multiprocessing.Pool(processes=4)                    # The multiprocessing pool initialization
            list1 = (vDummy[rbin], thermal, *fit_range)                 # The lists of fuctions and parameters are
            list2 = (vDummy[rbin], fit_dexp, *fit_range_dexp)           #   initialized for the each fit
            list3 = (vDummy[rbin], blastwave_inv, *fit_range)
            list4 = (vDummy[rbin], bosein, *fit_range)
            result = pool.starmap(apply_fit, [list1,list2,list3,list4]) # All fits are performed simultaneously in parallel
            pool.close()
            pool.join()

            out_thermal       = result[0] # Thermal
            out_fit_dexp      = result[1] # Double-exponential
            out_blastwave_inv = result[2] # Blast-Wave
            out_bosein        = result[3] # Bose-Einstein

            # Make the fit results fancy
            make_fancy_histogram(out_thermal,       25,  6, 3,  6, 0.5, 3)
            make_fancy_histogram(out_fit_dexp,      24,  4, 3,  4, 0.5, 3)
            make_fancy_histogram(out_blastwave_inv, 23,  1, 2, 14, 0.5, 3)
            make_fancy_histogram(out_bosein,        26, 42, 3, 42, 0.5, 3)

            # Save fit to the output file
            save_fit(out_thermal,       pname,             f'thermal_{common_suffix}_{rapidity_postfix}')
            save_fit(out_fit_dexp,      pname,  f'double_exponential_{common_suffix}_{rapidity_postfix}')
            save_fit(out_blastwave_inv, pname, f'blastwave_invariant_{common_suffix}_{rapidity_postfix}')
            save_fit(out_bosein,        pname,              f'bosein_{common_suffix}_{rapidity_postfix}')

            # Draw the fit result
            out_thermal.Draw('same PL')
            out_fit_dexp.Draw('same PL')
            out_blastwave_inv.Draw('same PL')
            out_bosein.Draw('same PL')

            x, y = array('d'), array('d')
            n = 0
            for hbin in range(1, vMC[rbin].GetSize() - 1):      # Here the BW fit is evaluated at the same points
                if vMC[rbin].GetBinContent(hbin):               #   as the Monte-Carlo histogram
                    n += 1                                      #   and then the TGraph with values
                    x.append(vMC[rbin].GetBinCenter(hbin))      #   y = Evaluated / Monte-Carlo is collected (vvv next line vvv)
                    y.append(out_blastwave_inv.Eval(vMC[rbin].GetBinCenter(hbin)) / vMC[rbin].GetBinContent(hbin))
            if n == 0:                                                # There is no data
                write_status(f'{function_name} {CURRENT_ST[order]}')  # Write the status and
                return                                                #   skip!
            gr = ROOT.TGraph(n, x, y)
            make_fancy_histogram(gr, 23,  1, 2, 14, 0.5, 3)
            ratio_cor.GetLowerPad().cd()                        # The ratio BW / Monte-Carlo
            gr.Draw('same LP')                                  # is drawn on the lower pad
            ratio_cor.GetUpperPad().cd()

            dndy_exp     = dndy_from_pt(vCorr[rbin], *dndy_exp_limits)             # dN/dy value from corrected spectra within defined pT range
            dndy_mc_low  = dndy_from_pt(vMC[rbin],   pt_start, dndy_exp_limits[0]) # dN/dy value from MC -- low pT part
            dndy_mc_high = dndy_from_pt(vMC[rbin],   dndy_exp_limits[1], pt_end)   # dN/dy value from MC -- high pT part

            blastwave.SetParameters(out_blastwave_inv.GetParameters())                     # BW fit for the d^2N/dptdy spectra
            save_fit(blastwave, pname, f'blastwave_pt_{common_suffix}_{rapidity_postfix}') # Fit parameters are taken from the 
                                                                                           # BW fit for the invariant pT spectra (above)
            tbox = ROOT.TPaveText(0.15, 0.05, 0.45, 0.5, 'NDC')                     # Text box with different information: fits, dN/dy etc
            bw_integral = blastwave.Integral(pt_start, pt_end)                      # Integral of the BW fit
            bw_high = blastwave.Integral(dndy_exp_limits[1], pt_end)                # Low pT part of the BW fit
            bw_low = blastwave.Integral(pt_start, dndy_exp_limits[0])               # High pT part of the BW fit
            tbox.AddText('BW fit range = {:.1f} - {:.1f} GeV/c'.format(*fit_range))
            tbox.AddText('BW Integral (0 - 6  GeV/c) = {:.3f}'.format(bw_integral))
            tbox.AddText('BW low p_{{T}}  = {:.3f} ({:.2f}%)'.format(bw_low, bw_low/bw_integral*100))
            tbox.AddText('BW high p_{{T}} = {:.3f} ({:.2f}%)'.format(bw_high, bw_high/bw_integral*100))
            tbox.AddText('dN/dy_{{exp}} [{:.2f}..{:.2f}] = {:.3e} \\pm {:.3e}'.format(*dndy_exp_limits, *dndy_exp))
            dndy_bw_low      = blastwave.Integral(pt_start, dndy_exp_limits[0])     # dN/dy value from BW fit -- low pT part
            dndy_bw_low_err  = math.fabs(dndy_bw_low - dndy_mc_low[0])              # dN/dy error from BW fit -- low pT part
            dndy_bw_high     = blastwave.Integral(dndy_exp_limits[1], pt_end)       # dN/dy value from BW fit -- high pT part
            dndy_bw_high_err = math.fabs(dndy_bw_high - dndy_mc_high[0])            # dN/dy error from BW fit -- high pT part
            tbox.AddText('dN/dy_{{BW}} [{:.2f}..{:.2f}] = {:.3e} \\pm {:.3e}'.format(0, dndy_exp_limits[0], dndy_bw_low , dndy_bw_low_err))
            tbox.AddText('dN/dy_{{BW}} [{:.2f}..{:.2f}] = {:.3e} \\pm {:.3e}'.format(dndy_exp_limits[1], pt_end, dndy_bw_high, dndy_bw_high_err))
            dndy_final     = dndy_bw_low + dndy_exp[0] + dndy_bw_high               # Final dN/dy value -- sum of previous dN/dy components
            dndy_final_err = math.sqrt(dndy_bw_low_err**2 + dndy_exp[1]**2 + dndy_bw_high_err**2) # final dN/dy error
            tbox.AddText('dN/dy_{{Final}} = {:.3e} \\pm {:.3e}'.format(dndy_final, dndy_final_err))
            # ~ tbox.Draw()
            particle_info = {}                   # A set of particle information (Python dict)
            particle_info['name'] = pname
            particle_info['centrality'] = cname
            particle_info['y'] = f'{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
            particle_info['bw_range'] = f'{fit_range[0]} - {fit_range[1]}'
            particle_info['bw'] = f'Integral: {bw_integral}, low pT: {bw_low}, high pT: {bw_high}'
            particle_info['dndy_bw_low'] = f'{dndy_bw_low} +/- {dndy_bw_low_err}  in  [{pt_start} - {dndy_exp_limits[0]}]'
            particle_info['dndy_exp'] = f'{dndy_exp[0]} +/- {dndy_exp[1]}  in  [{dndy_exp_limits[0]} - {dndy_exp_limits[1]}]'
            particle_info['dndy_bw_high'] = f'{dndy_bw_high} +/- {dndy_bw_high_err}  in  [{dndy_exp_limits[1]} - {pt_end}]'
            particle_info['dndy'] = f'{dndy_final} +/- {dndy_final_err}'
            if TBL_FILE is not None:
                with open(TBL_FILE, 'a') as f:                   # All information above will be written in the 
                    f.write(json.dumps(particle_info, indent=4)) #   output JSON file (TBL_FILE)
        legend.Draw()
        ratio_cor.GetLowerPad().cd()
        # ~ g_uncor.Draw()
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the \f$ dN/dy \f$ spectra
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param canvas TCanvas prepared for drawing
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
def draw_dndy(do_dedx, canvas, p_pos, cbin, fname):
    order = 2
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]                     # Particle name
    cbins    = list(Centrality)[cbin]                     # Centrality bin limits
    cname    = f'{cbins[0]} - {cbins[1]}%'                # Centrality bin name
    pt_start = list(Particles.values())[0]['pt_bins'][1]  # pT limits: low and 
    pt_end   = list(Particles.values())[0]['pt_bins'][2]  #   high
    print(f' draw dN/dy =============== {pname},  {cname},  {pt_start} < pT < {pt_end}  =============== ')

    if(do_dedx):
        postfix='_dedx'
    else:
        postfix=''
    canvas.cd()
    vMC   = []
    vCorr = []
    vBW   = []
    common_suffix = f'{pname}_centrality{cbin}'
    for rbin in range(len(rapidity_bins)):
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        common_postfix = f'{postfix}_{rapidity_postfix}'
        mc   = OUT_FILE.Get(f'{pname}/h__pt_{common_suffix}_mc_{rapidity_postfix}').Clone('mc')       # Monte-Carlo histogram
        corr = OUT_FILE.Get(f'{pname}/h__pt_{common_suffix}_corr{common_postfix}').Clone('corr')      # Reconstructed (corrected) histogram
        try:
            fit  = OUT_FILE.Get(f'{pname}/blastwave_pt_{common_suffix}_{rapidity_postfix}').Clone('fit')  # The Blast-Wave fit (not the invariant one!)
        except ReferenceError:
            print(f'  dN/dy:   no Blast-Wave fit found -- skipping')
            write_status(f'{function_name} {CURRENT_ST[order]}')
            return
        vMC.append(mc)
        vCorr.append(corr)
        vBW.append(fit)

    # The the integral of the reconstructed points within 
    #   the 'dndy_exp_limits = [low, high]' range will be taken 
    #   for the total dN/dy value calculation
    if (pname == 'pim' or pname == 'pip'):  # pi- and pi+
        dndy_exp_limits = [0.1, 1.5]
    elif (pname == 'km'):                   # K-
        dndy_exp_limits = [0.2, 1.5]
    elif (pname == 'kp'):                   # K+
        dndy_exp_limits = [0.2, 1.5]
    elif (pname == 'p'):                    # p
        dndy_exp_limits = [0.5, 1.5]
    elif (pname == 'd'):                    # d
        dndy_exp_limits = [0.8, 2.6]
    elif (pname == 't'):                    # t
        dndy_exp_limits = [1.2, 3.9]
    elif (pname == 'He3'):                  # 3He
        dndy_exp_limits = [1.2, 3.9]

    x_mc, y_mc   = array('d'), array('d')   # Vectors of the Monte-Carlo dN/dy spectra points values
    ex_mc, ey_mc = array('d'), array('d')   #   end errors
    x, y   = array('d'), array('d')         # Vectors of the reconstructed dN/dy spectra points values
    ex, ey = array('d'), array('d')         #   end errors
    for rbin in range(len(rapidity_bins)):
        dndy_exp         = dndy_from_pt(vCorr[rbin], *dndy_exp_limits)              # The recomstructed points part of the reconstructed dN/dy

        dndy_mc_low      = dndy_from_pt(vMC[rbin],   pt_start, dndy_exp_limits[0])  # The low-pT part of the Monte-Carlo dN/dy
        dndy_bw_low      = vBW[rbin].Integral(pt_start, dndy_exp_limits[0])         # The low-pT Blast-Wave fit part of the reconstructed dN/dy
        dndy_bw_low_err  = math.fabs(dndy_bw_low - dndy_mc_low[0])                  # The Blast-Wave fit error if the difference between the
                                                                                    #   BW fit value and the MC value within the same low-pT range
        dndy_mc_high     = dndy_from_pt(vMC[rbin],   dndy_exp_limits[1],   pt_end)  # The high-pT part of the Monte-Carlo dN/dy
        dndy_bw_high     = vBW[rbin].Integral(dndy_exp_limits[1],   pt_end)         # The high-pT Blast-Wave fit part of the reconstructed dN/dy
        dndy_bw_high_err = math.fabs(dndy_bw_high - dndy_mc_high[0])                # The Blast-Wave fit error if the difference between the
                                                                                    #   BW fit value and the MC value within the same high-pT range
        dndy_final       = dndy_bw_low + dndy_exp[0] + dndy_bw_high                 # The final reconstructed dN/dy = sum of the reconstructed points and
                                                                                    #   the Blast-Wave fit (low- and high-pT parts)
        dndy_final_err   = math.sqrt(dndy_bw_low_err**2 + dndy_exp[1]**2 + dndy_bw_high_err**2) # The final dN/dy error is calculated as the square root of the sum
                                                                                                #   of the squares of the BW fit errors and the reconstructed points integral error  
        x.append((rapidity_bins[rbin][0] + rapidity_bins[rbin][1]) / 2)
        y.append(dndy_final)
        ex.append(0)
        ey.append(dndy_final_err)
        
        dndy_orig      = dndy_from_pt(vMC[rbin],   pt_start, pt_end)
        x_mc.append((rapidity_bins[rbin][0] + rapidity_bins[rbin][1]) / 2)
        y_mc.append(dndy_orig[0])
        ex_mc.append(0)
        ey_mc.append(dndy_orig[1])
    n = len(x)
    gr = ROOT.TGraphErrors(n, x, y, ex, ey)
    make_fancy_histogram(gr, 20, 1, 1, 1, msize = 2)
    hmax = ROOT.TMath.MaxElement(n, gr.GetY()) * 1.5 # add 50% on the top
    hmin = ROOT.TMath.MinElement(n, gr.GetY()) * 0.5 # add 50% on the bottom
    if hmin < 0.001: hmin = 0
    frame = canvas.DrawFrame(-1.19, hmin, 1.19, hmax)
    make_fancy_frame(frame, f'{SYS}, {pname}, {cname}', 'y', 'dN/dy')
    gr.Draw('same PE')

    gr_mc = ROOT.TGraphErrors(n, x_mc, y_mc, ex_mc, ey_mc)
    make_fancy_histogram(gr_mc, 24, 2, 1, 2, msize = 2)
    gr_mc.Draw('same PE')

    canvas.Print(fname)
    canvas.Clear()
    canvas.SetLogy(0)

    grsrc = gr.Clone(f'h__rapidity_{common_suffix}')
    grsmc = gr_mc.Clone(f'h__rapidity_mc_{common_suffix}')
    OUT_FILE.cd(f'{pname}')
    grsrc.Write()
    grsmc.Write()
    OUT_FILE.cd()
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the coalescence parameter B2
#  \param canvas TCanvas prepared for drawing
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
def draw_b2(canvas, cbin, fname):
    order = 4
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    canvas.cd()
    rc_p = []  # Array of the proton pT spectra within different rapidity intervals
    rc_d = []  # Array of the deuteron pT spectra within different rapidity intervals

    for rbin in range(len(rapidity_bins)):
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hp   = OUT_FILE.Get(f'p/h__pt_p_centrality{cbin}_corr_{rapidity_postfix}').Clone('hp')
        hd   = OUT_FILE.Get(f'd/h__pt_d_centrality{cbin}_corr_{rapidity_postfix}').Clone('hd')
        rc_p.append(hp)
        rc_d.append(hd)

    p_exp_limits = [0.5, 1.5]                                 # The proton points will be taken within this range
    p_pt_start = list(Particles.values())[4]['pt_bins'][1]    # The original proton spectra low limit from the settings file
    p_pt_end   = list(Particles.values())[4]['pt_bins'][2]    # The original proton spectra high limit from the settings file
    p_pt_bins  = list(Particles.values())[4]['pt_bins'][0]    # The number of bins for the original proton spectra
    p_binw     = (p_pt_end - p_pt_start) / p_pt_bins          # The proton spectra bin width
    p_exp_bins = (p_exp_limits[1] - p_exp_limits[0]) / p_binw # The number of the proton spectra bins within the chosen range for the B2 analysis

    d_exp_limits = [0.8, 2.6]                                 # The deuteron points will be taken within this range
    d_pt_start = list(Particles.values())[6]['pt_bins'][1]    # The original deuteron spectra low limit from the settings file
    d_pt_end   = list(Particles.values())[6]['pt_bins'][2]    # The original deuteron spectra high limit from the settings file
    d_pt_bins  = list(Particles.values())[6]['pt_bins'][0]    # The number of bins for the original deuteron spectra
    d_binw     = (d_pt_end - d_pt_start) / d_pt_bins          # The deuteron spectra bin width -- !!! must be 2 times wider than the proton spectra bin width !!!
    d_exp_bins = (d_exp_limits[1] - d_exp_limits[0]) / d_binw # The number of the deuteron spectra bins within the chosen range for the B2 analysis

    b2_limits = [0, 0]                                                        # The coalescence parameter B2 analysis range
    if p_exp_limits[0]*2. < d_exp_limits[0]:  b2_limits[0] = d_exp_limits[0]  # Low limit: = 'deuteron low limit' if 'd low limit' > 'p low limit X2'
    else: b2_limits[0] = p_exp_limits[0]*2.                                   #            = 'proton low limit X2' otherwise
    if p_exp_limits[1]*2. > d_exp_limits[1]:  b2_limits[1] = d_exp_limits[1]  # High limit: = 'deuteron high limit' if 'd high limit' < 'p high limit X2'
    else: b2_limits[1] = p_exp_limits[1]*2.                                   #             = 'proton high limit X2' otherwise

    b2_bins = int((b2_limits[1] - b2_limits[0]) / d_binw)                     # The number of B2 spectra bins with the bin width of deuteron (= proton bins width X2)
    scale_f = 1E3                                                             # The scale factor for the spectra. Not really needed, but the plot Y axis labels became nicer
    for rbin in range(len(rapidity_bins)):                                    # For the pT spectra in the each rapidity bin:
        x, y   = array('d'), array('d')                                       #   Initialize B2 points values
        ex, ey = array('d'), array('d')                                       #   and errors arrays
        for b2bin in range(b2_bins):                                          #   Then, for the each B2 bin (point):
            d_llim = b2_limits[0] + d_binw*b2bin                              #     Calculate the low limit for the corresponding deuteron spectra bin (point)
            d_rlim = b2_limits[0] + d_binw + d_binw*b2bin                     #     Calculate the high limit for the corresponding deuteron spectra bin (point)
            p_llim = b2_limits[0]/2. + p_binw*b2bin                           #     Calculate the low limit for the corresponding proton spectra bin (point)
            p_rlim = b2_limits[0]/2. + p_binw + p_binw*b2bin                  #     Calculate the high limit for the corresponding proton spectra bin (point)
            nd = get_ptbin_content(rc_d[rbin], d_llim, d_rlim)                #     Get the value of in the selected deuteron spectra bin
            np = get_ptbin_content(rc_p[rbin], p_llim, p_rlim)                #     Get the value of in the selected proton spectra bin
            x.append(d_llim + d_binw/2.)                                      #     Set the deuteron pT bin center as the B2 point X value
            b2 = nd[0] / (np[0]*np[0]) * (2*3.14) * scale_f                   #     Calculate the B2 coalescence parameter applying the scale factor
            y.append(b2)                                                      #     Set the calculated parameter  as the B2 point Y value
            if nd[0]:                                                         #     Check if the value in the deuteron spectra bin (nd[0])  is not equal to 0
                d_rel_error = nd[1] / nd[0]                                   #     Calculate the relative error for the d bin (nd[1] = absolute error)
            else:
                d_rel_error = 0
            if np[0]:                                                         #     Check if the value in the proton spectra (np[0]) bin is not equal to 0
                p_rel_error = np[1] / np[0]                                   #     Calculate the relative error for the p bin (np[1] = absolute error)
            else:
                p_rel_error = 0
            p_mult_err  = p_rel_error + p_rel_error                           #     Multiplication error = sum of the relative errors
            b2_rel_err  = d_rel_error + p_mult_err                            #     Division error = sum of the relative errors
            b2_abs_err  = b2_rel_err * b2                                     #     Calculate the absolute error
            ex.append(0)
            ey.append(b2_abs_err)                                             #     Set the absolute error as the B2 point Y error
        n = len(x)
        gr = ROOT.TGraphErrors(n, x, y, ex, ey)                          #   Create the TGraphErrors object using the calculated points
        make_fancy_histogram(gr, 20, 1, 1, 1, msize = 2)
        hmax = ROOT.TMath.MaxElement(n, gr.GetY()) * 1.2 # add 20% on the top
        hmin = ROOT.TMath.MinElement(n, gr.GetY()) * 0.8 # add 20% on the bottom
        frame = canvas.DrawFrame(b2_limits[0]*0.9, hmin, b2_limits[1]*1.1, hmax)
        power = math.floor(math.log10(scale_f))
        make_fancy_frame(frame, f'{SYS}, {rapidity_bins[rbin][0]} < y < {rapidity_bins[rbin][1]}, {cname}', 
                        'p_{T}, GeV/c', f'B2 #times 10^{{{power}}}, GeV^{{2}}/c^{{3}}')
        gr.Draw('same PE')
        grb2 = gr.Clone(f'h__b2_centrality{cbin}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}')
        OUT_FILE.cd(f'd/')
        grb2.Write()
        OUT_FILE.cd()
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the coalescence parameter B3 for the He3
#  \param canvas TCanvas prepared for drawing
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
def draw_b3he3(canvas, cbin, fname):
    order = 10
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    canvas.cd()
    rc_p = []   # Array of the proton pT spectra within different rapidity intervals
    rc_he3 = [] # Array of the He3 pT spectra within different rapidity intervals

    for rbin in range(len(rapidity_bins)):
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hp   = OUT_FILE.Get(f'p/h__pt_p_centrality{cbin}_corr_{rapidity_postfix}').Clone('hp')
        try:
            hhe3 = OUT_FILE.Get(f'He3/h__pt_He3_centrality{cbin}_corr_{rapidity_postfix}').Clone('hhe3')
        except ReferenceError:
            print(f'  B3(He3):   no pT spectra found -- skipping')
            write_status(f'{function_name} {CURRENT_ST[order]}')
            return
        rc_p.append(hp)
        rc_he3.append(hhe3)

    p_exp_limits = [0.4, 1.3]                                 # Proton points will be taken within this range
    p_pt_start = list(Particles.values())[4]['pt_bins'][1]    # The original proton spectra low limit from the settings file
    p_pt_end   = list(Particles.values())[4]['pt_bins'][2]    # The original proton spectra high limit from the settings file
    p_pt_bins  = list(Particles.values())[4]['pt_bins'][0]    # The number of bins for the original proton spectra
    p_binw     = (p_pt_end - p_pt_start) / p_pt_bins          # The proton spectra bin width
    p_exp_bins = (p_exp_limits[1] - p_exp_limits[0]) / p_binw # The number of the proton spectra bins within the chosen range for the B2 analysis

    he3_exp_limits = [0.9, 3.9]                                       # The He3 points will be taken within this range
    he3_pt_start = list(Particles.values())[8]['pt_bins'][1]          # The original He3 spectra low limit from the settings file
    he3_pt_end   = list(Particles.values())[8]['pt_bins'][2]          # The original He3 spectra high limit from the settings file
    he3_pt_bins  = list(Particles.values())[8]['pt_bins'][0]          # The number of bins for the original He3 spectra
    he3_binw     = (he3_pt_end - he3_pt_start) / he3_pt_bins          # The He3 spectra bin width -- !!! must be 3 times wider than the proton spectra bin width !!!
    he3_exp_bins = (he3_exp_limits[1] - he3_exp_limits[0]) / he3_binw # The number of the He3 spectra bins within the chosen range for the B3 analysis

    b3_limits = [0, 0]                                                           # The coalescence parameter B3 analysis range
    if p_exp_limits[0]*3. < he3_exp_limits[0]:  b3_limits[0] = he3_exp_limits[0] # Low limit: = 'He3 low limit' if 'He3 low limit' > 'p low limit X3'
    else: b3_limits[0] = p_exp_limits[0]*3.                                      #            = 'proton low limit X3' otherwise
    if p_exp_limits[1]*3. > he3_exp_limits[1]:  b3_limits[1] = he3_exp_limits[1] # High limit: = 'He3 high limit' if 'He3 high limit' < 'p high limit X3'
    else: b3_limits[1] = p_exp_limits[1]*3.                                      #             = 'proton high limit X3' otherwise

    b3_bins = int((b3_limits[1] - b3_limits[0]) / he3_binw)                   # The number of B3 spectra bins with the bin width of He3 (= proton bins width X3)
    scale_f = 1E3                                                             # The scale factor for the spectra. Not really needed, but the plot Y axis labels became nicer
    for rbin in range(len(rapidity_bins)):                                    # For the pT spectra in the each rapidity bin:
        x, y   = array('d'), array('d')                                       #   Initialize B2 points values
        ex, ey = array('d'), array('d')                                       #   and errors arrays
        for b3bin in range(b3_bins):                                          #   Then, for the each B3 bin (point):
            he3_llim = b3_limits[0] + he3_binw*b3bin                          #     Calculate the low limit for the corresponding He3 spectra bin (point)
            he3_rlim = b3_limits[0] + he3_binw + he3_binw*b3bin               #     Calculate the high limit for the corresponding He3 spectra bin (point)
            p_llim   = b3_limits[0]/3. + p_binw*b3bin                         #     Calculate the low limit for the corresponding proton spectra bin (point)
            p_rlim   = b3_limits[0]/3. + p_binw + p_binw*b3bin                #     Calculate the high limit for the corresponding proton spectra bin (point)
            nhe3 = get_ptbin_content(rc_he3[rbin], he3_llim, he3_rlim)        #     Get the value of in the selected He3 spectra bin
            np = get_ptbin_content(rc_p[rbin], p_llim, p_rlim)                #     Get the value of in the selected proton spectra bin
            x.append(he3_llim + he3_binw/2.)                                  #     Set the He3 pT bin center as the B3 point X value
            b3 = nhe3[0] / (np[0]*np[0]*np[0]) * pow(2*3.14, 2) * scale_f     #     Calculate the B3(He3) coalescence parameter applying the scale factor
            y.append(b3)                                                      #     Set the calculated parameter as the B3 point Y value
            if nhe3[0]:                                                       #     Check if the value in the He3 spectra bin (nhe3[0]) is not equal to 0
                he3_rel_error = nhe3[1] / nhe3[0]                             #     Calculate the relative error for the He3 bin (nhe3[1] = absolute error)
            else:
                he3_rel_error = 0
            if np[0]:                                                         #     Check if the value in the proton spectra (np[0]) bin is not equal to 0
                p_rel_error = np[1] / np[0]                                   #     Calculate the relative error for the p bin (np[1] = absolute error)
            else:
                p_rel_error = 0
            p_mult_err  = p_rel_error + p_rel_error + p_rel_error             #     Multiplication error = sum of the relative errors
            b3_rel_err  = he3_rel_error + p_mult_err                          #     Division error = sum of the relative errors
            b3_abs_err  = b3_rel_err * b3                                     #     Calculate the absolute error
            ex.append(0)
            ey.append(b3_abs_err)                                             #     Set the absolute error as the B3 point Y error
        n = len(x)
        gr = ROOT.TGraphErrors(n, x, y, ex, ey)                         #   Create the TGraphErrors object using the calculated points
        make_fancy_histogram(gr, 20, 1, 1, 1, msize = 2)
        hmax = ROOT.TMath.MaxElement(n, gr.GetY()) * 1.2 # add 20% on the top
        hmin = ROOT.TMath.MinElement(n, gr.GetY()) * 0.8 # add 20% on the bottom
        frame = canvas.DrawFrame(b3_limits[0]*0.9, hmin, b3_limits[1]*1.1, hmax)
        power = math.floor(math.log10(scale_f))
        make_fancy_frame(frame, f'{SYS}, {rapidity_bins[rbin][0]} < y < {rapidity_bins[rbin][1]}, {cname}', 
                        'p_{T}, GeV/c', f'B3(^{{3}}He) #times 10^{{{power}}}, GeV^{{2}}/c^{{3}}')
        gr.Draw('same PE')
        grb3 = gr.Clone(f'h__b3_centrality{cbin}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}')
        OUT_FILE.cd(f'He3/')
        grb3.Write()
        OUT_FILE.cd()
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine makes the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra for the defined rapidity bins
#  \param do_corr Switch for the efficiency and contamination corrections:
#     - 1 = apply
#     - 0 = do not apply
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param particle The particle name (p, d, He4 etc)
#  \param cbin The centrality bin (0, 1, 2, etc)
#  \param psname The name of the selected phase space histogram in the input file
#  \param postfix The phase space histogram postfix, can be:
#     - 'mc' - For the Monte-Carlo phase space.
#     - 'uncor' - For the uncorrected reconstructed phase space.
#     - 'corr' - For the corrected reconstructed phase space.
def make_pt_spectra_mpdpid(do_corr, do_dedx, particle, cbin, psname, postfix):
    try:
        order = 0
        CURRENT_ST[order] += 1
        function_name = sys._getframe().f_code.co_name
        if skip_function(function_name, order): return

        cbins    = list(Centrality)[cbin]
        cname    = f'{cbins[0]} - {cbins[1]}%'

        if((postfix == 'uncorr' or postfix == 'corr') and do_dedx):
            postfix=f'{postfix}_dedx'
        elif ((postfix == 'uncorr' or postfix == 'corr') and not do_dedx):
            postfix=f'{postfix}'
        hPhaseSpace = INP_FILE.Get(psname).Clone(f'phasespace_{postfix}_{particle}_centrality{cbin}') # The selected phase-space
        hEvents = INP_FILE.Get('h__events').Clone('hEvents')                                          # The histogram with the events number for the each centrality bin
        name_ptspectra = f'h__pt_{particle}_centrality{cbin}_{postfix}'                               # This is the histogram name in the output ROOT file

        # To obtain the final reconstructed spectra the uncorrected phase-space must be corrected by the efficiencies and contamination
        if(do_corr): 
        # --- TPC, Secondaries
            # ------- efficiency
            apply_efficiency(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__efficiency_tpc_{particle}_centrality{cbin}').Clone('eff'))
            # ------- contamination
            apply_contamination(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__contamination_secondaries_{particle}_centrality{cbin}').Clone('cont'))

            if(do_dedx):
            # --- PID TPC only (dE/dx)
                # ------- efficiency
                apply_efficiency(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__efficiency_pid_dedx_{particle}_centrality{cbin}').Clone('eff'))
                # ------- contamination
                apply_contamination(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__contamination_pid_dedx_{particle}_centrality{cbin}').Clone('cont'))
            else:
            # --- ToF
                # ------- efficiency
                apply_efficiency(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__efficiency_tof_{particle}_centrality{cbin}').Clone('eff'))
            # --- PID combined
                # ------- efficiency
                apply_efficiency(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__efficiency_pid_{particle}_centrality{cbin}').Clone('eff'))
                # ------- contamination
                apply_contamination(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__contamination_pid_{particle}_centrality{cbin}').Clone('cont'))


        dn = hEvents.GetBinContent(cbin+1)      # The number of events in the selected centrality bin
        if dn == 0: return
        hPhaseSpace.Scale(1./ dn)               # The phase-space is also scaled by the number of events
        OUT_FILE.cd(f'{particle}')
        hPhaseSpace.Write()                         # The corrected phase-space is written to the output ROOT file
        for rbin in range(len(rapidity_bins)):      # Now, for the each rapidity bin
            rb_low = rapidity_bins[rbin][0]         #   the lower and
            rb_high = rapidity_bins[rbin][1]        #   the upper edge is selected,
            rb_width = math.fabs(rb_high - rb_low)  #   the rapidity bin width is calculated

            for hbin in range(1, hPhaseSpace.GetNbinsX()):                                             # In this loop the number
                if(math.fabs(hPhaseSpace.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA):      #   of the lower phase-space bin
                    left_edge = hbin                                                                   #   along the X-axis (rapidity) is calculated
                    break
            right_edge = int(left_edge + rb_width / hPhaseSpace.GetXaxis().GetBinWidth(left_edge) - 1) # Then the upper bin number
            res_name = f'{name_ptspectra}_y{rb_low}_{rb_high}'
            hResult = hPhaseSpace.ProjectionY(res_name, left_edge, right_edge)                         # The Y projection within these bin range ([lower, upper]) is selected
            hResult.Scale(1./ rb_width)                                                                #   and scaled by the rapidity bin width (otherwise it will be a sum)
            for k in range(1, hResult.GetSize() - 1):                                                  # Now, for the each bin of the projection:
                content = hResult.GetBinContent(k)                                                     #   N is extracted as the bin content
                if(content <= 0):                                                                      #   If there is no content in the bin
                    continue                                                                           #     skip it
                pt_binw = hResult.GetBinWidth(k)                                                       #   dpT is extracted as the bin width
                pt_mean = hResult.GetBinCenter(k)                                                      #   pT is extracted as the bin center
                error   = hResult.GetBinError(k)                                                       #   The bin error is extracted too
                hResult.SetBinContent(k, content / (pt_mean * pt_binw))                                #   The new bin content is calculated: N / pt * dpt * dy (dy -- on the previous step)
                hResult.SetBinError(k, error / (pt_mean * pt_binw))                                    #   as well as the new bin error
            hResult.GetXaxis().SetTitle('p_{T}, GeV/c')                                                # The X-axis title
            hResult.GetYaxis().SetTitle('d^{2}n/p_{T}dp_{T}dy')                                        # The Y-axis title
            hResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')                  # The histogram title
            hResult.Write()                                                                            # The histogram is written to the output ROOT file
        OUT_FILE.cd()
        write_status(f'{function_name} {CURRENT_ST[order]}')
    except:
        return



##  Same as postprocess.make_pt_spectra_mpdpid, but for the 'evPID' method
#  \param pid_type Switch for the evPID particle identification type:
#     - tpc = for the TPC dE/dx
#     - tof = for the ToF  \f$ m^{2} \f$
#  \param do_corr Switch for the efficiency and contamination corrections:
#     - 1 = apply
#     - 0 = do not apply
#  \param particle The particle name (p, d, He4 etc)
#  \param cbin The centrality bin (0, 1, 2, etc)
#  \param psname The name of the selected phase space histogram in the input file
def make_pt_spectra_evpid(pid_type, do_corr, particle, cbin, psname):
    order = 8
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    postfix = 'uncorr'
    if do_corr:
        postfix = 'corr'
    hPhaseSpace = INP_FILE.Get(psname).Clone(f'phasespace_{pid_type}_{postfix}_{particle}_centrality{cbin}') # The selected phase-space
    hEvents = INP_FILE.Get('h__events').Clone('hEvents')                                                     # The histogram with the events number for the each centrality bin
    name_ptspectra = f'h__pt_{pid_type}_{particle}_centrality{cbin}_{postfix}'                               # This is the histogram name in the output ROOT file

    # To obtain the final reconstructed spectra the uncorrected phase-space must be corrected by the efficiencies and contamination
    if do_corr:
        # --- TPC, Secondaries
        # ------- efficiency
        apply_efficiency(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__efficiency_tpc_{particle}_centrality{cbin}').Clone('eff'))
        # ------- contamination
        apply_contamination(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__contamination_secondaries_{particle}_centrality{cbin}').Clone('cont'))
        if pid_type == 'tof':
            # --- ToF
            # ------- efficiency
            apply_efficiency(hPhaseSpace, OUT_FILE.Get(f'{particle}/h__efficiency_tof_{particle}_centrality{cbin}').Clone('eff'))
        # --- PID
        # ------- efficiency
        hPhaseSpace.Scale(1./0.9544) # hardcoded -- 2 Sigma

    dn = hEvents.GetBinContent(cbin+1)      # The number of events in the selected centrality bin
    hPhaseSpace.Scale(1./ dn)               # The phase-space is also scaled by the number of events
    OUT_FILE.cd(f'{particle}')
    hPhaseSpace.Write()                         # The corrected phase-space is written to the output ROOT file
    for rbin in range(len(rapidity_bins)):      # Now, for the each rapidity bin
        rb_low = rapidity_bins[rbin][0]         #   the lower and
        rb_high = rapidity_bins[rbin][1]        #   the upper edge is selected,
        rb_width = math.fabs(rb_high - rb_low)  #   the rapidity bin width is calculated

        for hbin in range(1, hPhaseSpace.GetNbinsX()):                                             # In this loop the number
            if(math.fabs(hPhaseSpace.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA):      #   of the lower phase-space bin
                left_edge = hbin                                                                   #   along the X-axis (rapidity) is calculated
                break
        right_edge = int(left_edge + rb_width / hPhaseSpace.GetXaxis().GetBinWidth(left_edge) - 1) # Then the upper bin number
        res_name = f'{name_ptspectra}_y{rb_low}_{rb_high}'
        hResult = hPhaseSpace.ProjectionY(res_name, left_edge, right_edge)                         # The Y projection within these bin range ([lower, upper]) is selected
        hResult.Scale(1./ rb_width)                                                                #   and scaled by the rapidity bin width (otherwise it will be a sum)
        for k in range(1, hResult.GetSize() - 1):                                                  # Now, for the each bin of the projection:
            content = hResult.GetBinContent(k)                                                     #   N is extracted as the bin content
            if(content <= 0):                                                                      #   If there is no content in the bin
                continue                                                                           #     skip it
            pt_binw = hResult.GetBinWidth(k)                                                       #   dpT is extracted as the bin width
            pt_mean = hResult.GetBinCenter(k)                                                      #   pT is extracted as the bin center
            error   = hResult.GetBinError(k)                                                       #   The bin error is extracted too
            hResult.SetBinContent(k, content / (pt_mean * pt_binw))                                #   The new bin content is calculated: N / pt * dpt * dy (dy -- on the previous step)
            hResult.SetBinError(k, error / (pt_mean * pt_binw))                                    #   as well as the new bin error
        hResult.GetXaxis().SetTitle('p_{T}, GeV/c')                                                # The X-axis title
        hResult.GetYaxis().SetTitle('d^{2}n/p_{T}dp_{T}dy')                                        # The Y-axis title
        hResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')                  # The histogram title
        hResult.Write()                                                                            # The histogram is written to the output ROOT file
    OUT_FILE.cd()
    write_status(f'{function_name} {CURRENT_ST[order]}')



## Similar to the postprocess.draw_pt_spectra_mpdpid, this subroutine draws 
#    the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra for the defined rapidity bins.
#  \param canvas TCanvas prepared for drawing
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
def draw_pt_spectra_evpid(canvas, p_pos, cbin, fname):
    order = 9
    function_name = sys._getframe().f_code.co_name
    CURRENT_ST[order] += 1
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]                         # Particle name
    cbins    = list(Centrality)[cbin]                         # Centrality bin limits
    cname    = f'{cbins[0]} - {cbins[1]}%'                    # Centrality bin name
    pt_start = list(Particles.values())[p_pos]['pt_bins'][1]  # pT limits: low and 
    pt_end   = list(Particles.values())[p_pos]['pt_bins'][2]  # high
    p_mass   = float(list(Particles.values())[p_pos]['Mass']) # Particle mass
    print(f' draw pT =============== {pname},  {cname},  {pt_start} < pT < {pt_end}  =============== ')

    canvas.cd()
    vMC = []                                                  # The list of Monte-Carlo histograms
    vUncorr_tpc = []                                          # The list of uncorrected histograms for the evPID TPC only PID
    vCorr_tpc = []                                            # The list of corrected histograms for the evPID TPC only PID
    vUncorr_tof = []                                          # The list of uncorrected histograms for the evPID ToF only PID
    vCorr_tof = []                                            # The list of corrected histograms for the evPID ToF only PID
    common_suffix = f'{pname}_centrality{cbin}'               # Histograms name suffix for each particle and centrality
    for rbin in range(len(rapidity_bins)):                    # Populate lists with histograms for the eah rapidity bin
        common_postfix = f'_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hMC     = OUT_FILE.Get(f'{pname}/h__pt_{common_suffix}_mc_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}').Clone('hMC')
        hUncorr_tpc = OUT_FILE.Get(f'{pname}/h__pt_tpc_{common_suffix}_uncorr{common_postfix}').Clone('hUncorr')
        hCorr_tpc   = OUT_FILE.Get(f'{pname}/h__pt_tpc_{common_suffix}_corr{common_postfix}').Clone('hCorr')
        hUncorr_tof = OUT_FILE.Get(f'{pname}/h__pt_tof_{common_suffix}_uncorr{common_postfix}').Clone('hUncorr')
        hCorr_tof   = OUT_FILE.Get(f'{pname}/h__pt_tof_{common_suffix}_corr{common_postfix}').Clone('hCorr')
        vMC.append(hMC)
        vUncorr_tpc.append(hUncorr_tpc)
        vCorr_tpc.append(hCorr_tpc)
        vUncorr_tof.append(hUncorr_tof)
        vCorr_tof.append(hCorr_tof)

    plot_range = [0, 6]                     # Common histograms plotting range
    if (pname == 'pim' or pname == 'pip'):  # pi- and pi+ parameters
        plot_range = [0, 1.59]              # Plot range for the selected particle
    elif (pname == 'km' or pname == 'kp'):  # K- and K+ -- same as above
        plot_range = [0, 1.69]
    elif (pname == 'p'):                    # p -- same as above
        plot_range = [0, 2.19]
    elif (pname == 'd'):                    # d -- same as above
        plot_range = [0, 3.09]
    elif (pname == 'He3'):                  # 3He -- same as above
        plot_range = [0, 3.99]

    canvas.SetLogy(1)                       # Set the canvas Y-axis as log

    legend = TLegend(0.65, 0.6, 0.9, 0.9)   # Create, modify and place legend
    legend.SetNColumns(1)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)

    for rbin in range(len(rapidity_bins)):
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'  # Output file name for the each plot
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'

        ratio_cor_tpc   = TRatioPlot(vCorr_tpc[rbin], vMC[rbin])                     # Corrected ratio plot
        ratio_uncor_tpc = TRatioPlot(vUncorr_tpc[rbin], vMC[rbin])                   # Uncorrected ratio plot
        ratio_cor_tof   = TRatioPlot(vCorr_tof[rbin], vMC[rbin])                     # Corrected ratio plot
        ratio_uncor_tof = TRatioPlot(vUncorr_tof[rbin], vMC[rbin])                   # Uncorrected ratio plot
        ratio_cor_tpc.SetSeparationMargin(0.0)                                       # The margin between upper and lower (ratio) part
        ratio_cor_tpc.SetRightMargin(0.01)                                           # The right margin
        ratio_cor_tpc.SetLeftMargin(0.12)                                            # The left margin

        make_fancy_histogram(vMC[rbin], MS_MC, LC_MC, 1, LC_MC)                      # Markers, lines and colours for histograms:
        make_fancy_histogram(vUncorr_tpc[rbin], 28, LC_UNCOR, 1, LC_UNCOR)     # MC, uncorrected and corrected
        make_fancy_histogram(vCorr_tpc[rbin], 34, LC_COR, 1, LC_COR)
        make_fancy_histogram(vUncorr_tof[rbin], 32, LC_UNCOR, 1, LC_UNCOR)
        make_fancy_histogram(vCorr_tof[rbin], MS_COR+1, LC_COR, 1, LC_COR)

        if(rbin == 0):                                                    # The legend must be created once
            legend.AddEntry(vMC[rbin],         'MC',                'P')  # but then it can be used for each
            legend.AddEntry(vUncorr_tpc[rbin], 'evPID TPC',        'PL')  # rapidity bin
            legend.AddEntry(vCorr_tpc[rbin],   'evPID TPC (corr)', 'PL')
            legend.AddEntry(vUncorr_tof[rbin], 'evPID ToF',        'PL')
            legend.AddEntry(vCorr_tof[rbin],   'evPID ToF (corr)', 'PL')

        ratio_uncor_tpc.Draw()                                        # Draw uncorrected plot to have this object in the future
        ratio_cor_tof.Draw()                                          # Draw uncorrected plot to have this object in the future
        ratio_uncor_tof.Draw()                                        # Draw uncorrected plot to have this object in the future
        g_uncor_tpc = ratio_uncor_tpc.GetLowerRefGraph()
        g_cor_tof   = ratio_cor_tof.GetLowerRefGraph()
        g_uncor_tof = ratio_uncor_tof.GetLowerRefGraph()
        make_fancy_histogram(g_uncor_tpc, 26, LC_UNCOR, 1, LC_UNCOR)  # evPID TPC uncorrected
        make_fancy_histogram(g_uncor_tof, 32, LC_UNCOR, 1, LC_UNCOR)  # evPID ToF uncorrected
        make_fancy_histogram(g_cor_tof, MS_COR+1, LC_COR, 1, LC_COR)  # evPID ToF corrected


        ratio_cor_tpc.Draw()                                        # Draw corrected evPID TPCresults, all futher histograms will be
        ROOT.gPad.Update()                                          #   plotted on these pads (lower and upper)
        ratio_cor_tpc.GetLowerRefGraph().SetLineColor(LC_COR)       # Set the corrected results line colour
        ratio_cor_tpc.GetLowerRefGraph().SetLineWidth(2)            #   and the line width
        make_fancy_histogram(ratio_cor_tpc.GetLowerRefGraph(), MS_COR, LC_COR, 1, LC_COR)

        y_low = ratio_cor_tpc.GetUpperPad().GetFrame().GetY1()      # The upper pad Y-axis maximum
        y_max = ratio_cor_tpc.GetUpperPad().GetFrame().GetY2()      #   and minimum
        ratio_cor_tpc.GetLowerRefGraph().SetMinimum(0)              # Set the lower pad minimum to 0
        ratio_cor_tpc.GetLowerRefGraph().SetMaximum(2)              # Set the lower par maximum to 2
        ratio_cor_tpc.GetLowerRefXaxis().SetRangeUser(*plot_range)  # Set the lower (and tus upper too) pad the defined above plot range
        if(p_pos == 7):                                  # Special case for tritons
            ratio_cor_tpc.GetLowerRefGraph().SetMaximum(4)
        if(p_pos >= 6):                                  # Special case for d
            ratio_cor_tpc.GetUpperRefYaxis().SetRangeUser(math.pow(10, y_max - 3), math.pow(10, y_max+0.2)) # 3 orders of magnitude less
        if(p_pos >= 8):                                  # Special case for He3 and He4
            ratio_cor_tpc.GetUpperRefYaxis().SetRangeUser(math.pow(10, y_max - 3.5), math.pow(10, y_max+0.5)) # 3.5 orders of magnitude less
        ratio_cor_tpc.GetLowerRefYaxis().SetNdivisions(6)
        ratio_cor_tpc.GetUpperPad().cd()                        # Go to the upper pad and draw histograms and fits there
        vMC[rbin].Draw('same hist PE')                          # Monte-Carlo histogram for the current rapidity bin
        vUncorr_tpc[rbin].Draw('same hist PE')                  # Uncorrected results histogram for the current rapidity bin
        vCorr_tpc[rbin].Draw('same hist PE')                    # Corrected results histogram for the current rapidity bin
        vUncorr_tof[rbin].Draw('same hist PE')                  # Uncorrected results histogram for the current rapidity bin
        vCorr_tof[rbin].Draw('same hist PE')                    # Corrected results histogram for the current rapidity bin

        legend.Draw()
        ratio_cor_tpc.GetLowerPad().cd()
        g_uncor_tpc.Draw('P')
        g_cor_tof.Draw('P')
        g_uncor_tof.Draw('P')
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine makes the overall efficiency (TPC + ToF + PID)
#     and contamination (secondaries + PID) histograms.
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param particle The particle name (p, d, He4 etc)
#  \param cbin The centrality bin (0, 1, 2, etc)
def make_overall_efficiency(do_dedx, particle, cbin):
    try:
        order = 5
        CURRENT_ST[order] += 1
        function_name = sys._getframe().f_code.co_name
        if skip_function(function_name, order): return
        cbins    = list(Centrality)[cbin]
        cname    = f'{cbins[0]} - {cbins[1]}%'

        # ------- efficiency
        eff = OUT_FILE.Get(f'{particle}/h__efficiency_tpc_{particle}_centrality{cbin}').Clone('EFF')  # The TPC efficiency
        overall_efficiency = eff.Clone(f'h__overall_efficiency_{particle}_centrality{cbin}')
        eff = OUT_FILE.Get(f'{particle}/h__efficiency_tof_{particle}_centrality{cbin}').Clone('EFF')  # The ToF efficiency
        overall_efficiency.Multiply(eff)
        eff = OUT_FILE.Get(f'{particle}/h__efficiency_pid_{particle}_centrality{cbin}').Clone('EFF')  # The PID efficiency
        overall_efficiency.Multiply(eff)

        # ------- contamination
        cont  = OUT_FILE.Get(f'{particle}/h__contamination_secondaries_{particle}_centrality{cbin}').Clone('CONT')  # The secondaries contamination
        overall_contamination = cont.Clone(f'h__overall_contamination_{particle}_centrality{cbin}')
        cont  = OUT_FILE.Get(f'{particle}/h__contamination_pid_{particle}_centrality{cbin}').Clone('CONT')          # The PID contamination
        overall_contamination.Multiply(cont)

        OUT_FILE.cd(f'{particle}')
        overall_efficiency.Write()     # The overall efficiency is written to the output ROOT file
        overall_contamination.Write()  # The overall contamination is written to the output ROOT file
        # The overall efficiency for the ach rapidity bin
        for rbin in range(len(rapidity_bins)):      # So, for the each rapidity bin
            rb_low = rapidity_bins[rbin][0]         #   the lower and
            rb_high = rapidity_bins[rbin][1]        #   the upper edge is selected,
            rb_width = math.fabs(rb_high - rb_low)  #   the rapidity bin width is calculated

            for hbin in range(1, overall_efficiency.GetNbinsX()):                                             # In this loop the number
                if(math.fabs(overall_efficiency.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA):      #   of the lower phase-space bin
                    left_edge = hbin                                                                          #   along the X-axis (rapidity) is calculated
                    break
            right_edge = int(left_edge + rb_width / overall_efficiency.GetXaxis().GetBinWidth(left_edge) - 1) # Then the upper bin number
            res_name_eff  = f'overall_efficiency_{particle}_centrality{cbin}_y{rb_low}_{rb_high}'
            res_name_cont = f'overall_contamination_{particle}_centrality{cbin}_y{rb_low}_{rb_high}'
            eResult = overall_efficiency.ProjectionY(res_name_eff, left_edge, right_edge)     # The Y projection within these bin range ([lower, upper]) is selected
            eResult.Scale(1. / (right_edge - left_edge))                                      #   and scaled by the number of projected bins
            eResult.GetXaxis().SetTitle('p_{T}, GeV/c')                                       # The projection X-axis title
            eResult.GetYaxis().SetTitle('Overall efficiency')                                 # The projection Y-axis title
            eResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')         # The projection title
            eResult.Write()                                                                   # The projection is written to the outpur ROOT file
            cResult = overall_contamination.ProjectionY(res_name_cont, left_edge, right_edge) # Same for the overall contamination histograms
            cResult.Scale(1. / (right_edge - left_edge))
            cResult.GetXaxis().SetTitle('p_{T}, GeV/c')
            cResult.GetYaxis().SetTitle('Overall contamination')
            cResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')
            cResult.Write()
        OUT_FILE.cd()
        write_status(f'{function_name} {CURRENT_ST[order]}')
    except:
        return



##  This subroutine draws the overall efficiency histograms
#  \param canvas TCanvas prepared for drawing
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
def draw_overall_efficiency(canvas, p_pos, cbin, fname):
    order = 6
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]
    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    canvas.cd()
    vEff = []
    for rbin in range(len(rapidity_bins)):
        common_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hEff     = OUT_FILE.Get(f'{pname}/overall_efficiency_{pname}_centrality{cbin}_{common_postfix}').Clone('hEff')
        vEff.append(hEff)

    plot_range = [0, 6]
    if (pname == 'pim' or pname == 'pip'):  # pi- and pi+
        plot_range = [0, 1.59]
    elif (pname == 'km' or pname == 'kp'):  # K- and K+
        plot_range = [0, 1.69]
    elif (pname == 'p'):                    # p
        plot_range = [0, 2.29]
    elif (pname == 'd'):                    # d
        plot_range = [0, 3.09]
    elif (pname == 'He3'):                  # 3He
        plot_range = [0, 4.29]

    for rbin in range(len(rapidity_bins)):
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
        make_fancy_histogram(vEff[rbin], MS_MC, LC_MC, 1, LC_MC, msize = 2)
        frame = canvas.DrawFrame(plot_range[0], 0, plot_range[1], 1.19)
        make_fancy_frame(frame, f'{SYS}, {pname}, {cname}, {rapidity_bins[rbin][0]} < y < {rapidity_bins[rbin][1]}',
                         'p_{T}, GeV/c', 'Overall efficiency')
        vEff[rbin].Draw('same hist P')
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the frame for the final histograms
#  \param frame The TH1F frame
#  \param title The title of the frame (histogram)
#  \param x_title The title of x-axis
#  \param y_title The title of y-axis
def make_fancy_frame(frame, title, x_title, y_title):
    frame.SetTitle(title)     # The frame title
    xax = frame.GetXaxis()    # The frame X-axis
    xax.SetTitle(x_title)     # The frame X-axis title
    xax.SetTitleFont(42)      # The frame X-axis title font
    xax.SetTitleSize(0.06)    # The frame X-axis title size
    xax.SetTitleOffset(1.05)  # The frame X-axis title offset
    xax.SetLabelFont(42)      # The frame X-axis label font
    xax.SetLabelSize(0.06)    # The frame X-axis label size
    xax.SetLabelOffset(0.01)  # The frame X-axis label offset
    yax = frame.GetYaxis()    # The frame Y-axis
    yax.SetTitle(y_title)     # The frame Y-axis title
    yax.SetTitleFont(42)      # The frame Y-axis title font
    yax.SetTitleSize(0.06)    # The frame Y-axis title size
    yax.SetTitleOffset(1.1)   # The frame Y-axis title offset
    yax.SetLabelFont(42)      # The frame Y-axis label font
    yax.SetLabelSize(0.06)    # The frame Y-axis label size
    yax.SetLabelOffset(0.01)  # The frame Y-axis label offset
    yax.SetNdivisions(8)      # The frame Y-axis number of divisions



##  This makes the histograms fancy
#  \param hist The histogram to modify
#  \param mstyle The marker style
#  \param mcolor The marker colour
#  \param lstyle The line style
#  \param lcolor The line colour
#
#  Optionally one can set
#  \param msize The size of the marker
#  \param lwidth The width of the line
def make_fancy_histogram(hist, mstyle, mcolor, lstyle, lcolor, msize = None, lwidth = None):
    hist.SetMarkerStyle(mstyle)   # Set the histogram marker style
    hist.SetMarkerColor(mcolor)   # Set the histogram marker colour
    if msize:                     # If the marker size 'msize' is defined
        hist.SetMarkerSize(msize) #   set the histogram marker size
    hist.SetLineStyle(lstyle)     # Set the histogram line style
    hist.SetLineColor(lcolor)     # Set the histogram line colour
    if lwidth:                    # If the line width 'lwidth' is defined
        hist.SetLineWidth(lwidth) #   set the histogram line width



##  This subroutine calculates the efficiency (contamination) correction
#  \param num The name of the selected efficiency (or contamination) numerator in the input file
#  \param denom The name of the selected efficiency (or contamination) denominator in the input file
#  \param res The resulting histogram name
#  \param particle The particle name (p, d, He4 etc)
def calculate_efficiency(num, denom, res, particle):
    try:
        hNumerator   = EFF_FILE.Get(num).Clone('hNumerator')      # The efficiency numerator
        hDenominator = EFF_FILE.Get(denom).Clone('hDenominator')  # The efficiency denominator
        hEfficiency  = hNumerator.Clone(res)                      # The efficiency histogram initialization
        hEfficiency.Divide(hDenominator)                          # The efficiency is calculated as 'Efficiency = Numerator / Denominator'
        OUT_FILE.cd(f'{particle}')                                # Go to the current particle folder in the output ROOT file
        hEfficiency.Write()                                       # Write the efficiency histogram
        OUT_FILE.cd()                                             # Go back to the '/' of the output ROOT file
    except:
        return



##  This subroutine applies the efficiency correction: phase-space / efficiency
#  \param phasespace The phase space to correct
#  \param eff The efficiency
def apply_efficiency(phasespace, eff):
    phasespace.Divide(eff) # the usual ROOT way



##  This subroutine applies the contamination correction as: phase-space * (1 - contamination)
#  \param phasespace The phase space to correct
#  \param cont The contamination
def apply_contamination(phasespace, cont):
    for i in range(1, cont.GetNbinsX()):                                      # For the each bin along the X-axis of the selected phase-space (rapidity)
        for j in range(1, cont.GetNbinsY()):                                  #   then for the each bin along the Y-axis of the selected phase-space (pT)
            cont_content = cont.GetBinContent(i, j)                           #     retrieve the bin(i,j) content of the contamination histogram
            ps_content   = phasespace.GetBinContent(i, j)                     #     and same for the phase-space histogram. Then replace the bin(i,j) content of the
            phasespace.SetBinContent(i, j, ps_content * (1. - cont_content))  #     phase-space histogram by the: old_phasespace_content * (1 - contamination_content)



##  This subroutine saves the fit function
#  \param fit The fit function object (e.g. TF1)
#  \param pname The particle name
#  \param name The name of the fit function
def save_fit(fit, pname, name):
    res = fit.Clone(name)
    OUT_FILE.cd(f'{pname}')
    res.Write()
    OUT_FILE.cd()



##  This subroutine calculates the dN/dy for the selected \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param hist The \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param pt_low The lower edge of the first \f$ p_{T} \f$ bin
#  \param pt_high The upper edge of the last \f$ p_{T} \f$ bin
#  \return The dN/dy value (dndy[0]) and error (dndy[1])
def dndy_from_pt(hist, pt_low, pt_high):
    dndy   = [0, 0]
    counts = 0
    err    = 0
    for h_bin in range(1, hist.GetNbinsX()+1):                                                      # In this loop
        if(math.fabs(hist.GetBinLowEdge(h_bin) - pt_low) < SMALL_DELTA):                              #   the number of the left bin
            left_edge = h_bin                                                                         #   is found
        if(math.fabs(hist.GetBinLowEdge(h_bin) + hist.GetBinWidth(h_bin)  - pt_high) < SMALL_DELTA):  #   and then the number of the right bin
            right_edge = h_bin                                                                        #   is found too
            break
    for h_bin in range(left_edge, right_edge+1):  # Then within that found bins range
        n   = hist.GetBinContent(h_bin)           #   the bin content N is taken
        pt  =  hist.GetBinCenter(h_bin)           #   the mean pT is taken as the bin center
        dpt =  hist.GetBinWidth(h_bin)            #   the dpT is taken as the bin width
        e   =  hist.GetBinError(h_bin)            #   the error is also taken
        if n:                                     #   If the bin content is not zero, 
            counts += n * pt * dpt                #   then add the dN/dy value (n * pt * dpt) to the dN/dy sum
            err += math.pow(e * pt, 2)            #   and add the squared error to the sum of the errors
    dndy[0] = counts                              # The final dN/dy value (dndy[0]) is the sum of the dN/dy values of all bins in the range
    dndy[1] = math.sqrt(err)                      # The final error (dndy[1]) is the squared root from the sum of squares
    return dndy 



##  This subroutine returns the selected bin content of the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param hist The \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param pt_low The lower edge of the \f$ p_{T} \f$ bin
#  \param pt_high The upper edge of the \f$ p_{T} \f$ bin
#  \return The \f$ p_{T} \f$ bin content (res[0]) and error (res[1])
def get_ptbin_content(hist, pt_low, pt_high):
    res = [0, 0]
    for h_bin in range(1, hist.GetNbinsX() + 1):                                                      # In this loop
        if(math.fabs(hist.GetBinLowEdge(h_bin) - pt_low) < SMALL_DELTA):                              #   the number of the left bin
            left_edge = h_bin                                                                         #   is found
        if(math.fabs(hist.GetBinLowEdge(h_bin) + hist.GetBinWidth(h_bin)  - pt_high) < SMALL_DELTA):  #   and then the number of the right bin
            right_edge = h_bin                                                                        #   is found too
            break
    for h_bin in range(left_edge, right_edge+1):  # Then within that found bin
        n   = hist.GetBinContent(h_bin)           #   the bin content N is taken
        e   = hist.GetBinError(h_bin)             #   the error is taken too
        if n:                                     #   and then, of the bin content is not equal to 0
            res[0] = n                            #   bin content and error are returned as the result
            res[1] = e
    return res



##  This subroutine converts the selected \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra to the \f$ d^{2}N/dp_{T}dy \f$
#  \param hist The \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra to convert
def convert_to_pt(hist):
    for h_bin in range(1, hist.GetNbinsX()+1):
        n   = hist.GetBinContent(h_bin)       # Get the bin content N
        pt  =  hist.GetBinCenter(h_bin)       # Get the bin mean pT
        e   =  hist.GetBinError(h_bin)        # Get the bin error
        if n:
            hist.SetBinContent(h_bin, n * pt) # Replace the bin content by the 'n * pT'
            hist.SetBinError(h_bin, e * pt)   # And recalculate the bin error



##  This subroutine checks if the file exists
#  \param fname The file name to check
def check_file_exist(fname):
    if(not os.path.isfile(fname)):
        print(f'--- {fname} does not exist')
        exit(-1)



##  This subroutine checks if the ROOT file can be opened and is not 'zombie'
#  \param file The ROOT file to check
def check_root_file(file):
    if (not file or file.IsZombie()):
        print(f'--- Can not open {file}')
        exit(-1)



##  This subroutine reads the saved progress of the funtion from the ASCII file (STATUS_FILE)
#  \param func The function name
#  \return The saved progress of the selected function
def read_status(func):
    res = ['Empty', 0]
    if os.path.exists(STATUS_FILE):                     # If the STATUS_FILE exists
        with open(STATUS_FILE, 'r') as f:               #   open it
            for line in f:
                if f'{func}' in line:                   #   find the line with the function name
                    res = line.rstrip('\n').split(' ')  #   and read the value (the function calls counter)
                    break
    else:
        with open(STATUS_FILE, 'w') as document: pass   # If file does not exists -- create it
    return(int(res[1]))



##  This subroutine checks if the function must be skipped.
#  \param fnc_name The function name
#  \param order The function execution order
#  \return True if the current function call counter is less or equal than the saved progress (skip),
#          false if not (do not skip)
def skip_function(fnc_name, order):
    saved_state = read_status(fnc_name)
    global CURRENT_ST
    global SKIP
    if SKIP and CURRENT_ST[order] <= saved_state:
        return True
    return False



##  This subroutine saves the current progress to the ASCII file (STATUS_FILE)
#  \param progress String of the format '(string)function_name (integer)counter' e.g.: ```draw_pt_something 2```
def write_status(progress):
    func, position = progress.rstrip('\n').split(' ')
    do_update = 0
    with open(STATUS_FILE, 'r+') as file:                                                     # Open the file
        for line in file:
            if f'{func}' in line:                                                             # If the function name is already here
                do_update = 1                                                                 #   skip to the 'update' part of this function
                break
        else:                                                                                 # Else
            file.write(f'{func} {position}\n')                                                #   write the function name and a calls counter
    if do_update:                                                                             # If do update
        for line in fileinput.input(STATUS_FILE, inplace=True):                               #   find the corresponding function line
            print(line.replace(f'{func} {int(position) - 1}', f'{func} {position}'), end='')  # and update the counter



##  This subroutine generates the pdf file with the report
#  \param is_extended Switch for the extended QA plots production:
#     - 1 = on
#     - 0 = off
#  \param dirname The input directory with produced plots to include into the report
#  \param fname The output report file (without the extention)
#
def generate_report(is_extended, dirname, fname):
    geometry_options = {'tmargin': '1cm', 'lmargin': '1cm', 'rmargin': '1cm', 'bmargin': '1.5cm'}
    doc = Document(geometry_options=geometry_options)
    doc.packages.append(Package('hyperref'))
    doc.packages.append(Package('graphicx'))
    cf.active = cf.Version1(indent=False)
    doc.append(Command('tableofcontents'))
    doc.append(NoEscape(r'\clearpage'))
    for part in Particles:
        with doc.create(Section(f'{part}')):
            if is_extended:
                list_efficiencies   = ['tpc', 'tof', 'pid']
                list_contaminations = ['secondaries', 'pid', 'tof']
                if DEDX:
                    list_efficiencies.append('pid_dedx')
                    list_contaminations.append('pid_dedx')
                for eff in list_efficiencies:
                    subsname=eff.upper()
                    if(eff == 'pid_dedx'):
                        subsname = 'PID dE/dx only'
                    with doc.create(Subsection(f'{subsname} efficiency')):
                        doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                        for cbin in range(0, len(Centrality)):
                            doc.append(Command('includegraphics', NoEscape(f'{dirname}/efficiency/{eff}/{part}_centrality{cbin}.{EXT}')))
                        doc.append(NoEscape(r'}'))
                doc.append(NoEscape(r'\newpage')) 
                for eff in list_contaminations:
                    subsname=eff.upper()
                    if(eff == 'pid_dedx'):
                        subsname = 'PID dE/dx only'
                    with doc.create(Subsection(f'{subsname} contamination')):
                        doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                        for cbin in range(0, len(Centrality)):
                            doc.append(Command('includegraphics', NoEscape(f'{dirname}/contamination/{eff}/{part}_centrality{cbin}.{EXT}')))
                        doc.append(NoEscape(r'}'))
                doc.append(NoEscape(r'\newpage')) 
            with doc.create(Subsection('Results')):
                doc.append(NoEscape(r'\begin{center}'))
                for cbin in range(0, len(Centrality)):
                    doc.append(NoEscape(r'\resizebox{0.94\textwidth}{!} {'))
                    for rbin in range(len(rapidity_bins)):
                        spectra = f'{dirname}/results/pt/{part}/{part}_centrality{cbin}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
                        doc.append(Command('includegraphics', NoEscape(spectra)))
                    doc.append(NoEscape(r'}\\'))
                doc.append(NoEscape(r'\end{center}'))
                if part != 'ap' and part != 'He4':
                    doc.append(NoEscape(r'\begin{center}'))
                    doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                    for cbin in range(0, len(Centrality)):
                        spectra = f'{dirname}/results/dndy/{part}/{part}_centrality{cbin}.{EXT}'
                        doc.append(Command('includegraphics', NoEscape(spectra)))
                    doc.append(NoEscape(r'}'))
                    doc.append(NoEscape(r'\end{center}'))
                if part == 'd':
                    doc.append(NoEscape(r'\begin{center}'))
                    for cbin in range(0, len(Centrality)):
                        doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                        for rbin in range(len(rapidity_bins)):
                            spectra = f'{dirname}/results/coal/b2_centrality{cbin}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
                            doc.append(Command('includegraphics', NoEscape(spectra)))
                        doc.append(NoEscape(r'}\\'))
                    doc.append(NoEscape(r'\end{center}'))
                doc.append(NoEscape(r'\newpage'))    
        doc.append(NoEscape(r'\clearpage'))
    doc.generate_pdf(fname, clean_tex=False)


# This is a truly main function
if __name__ == '__main__':
    main()
    print('All done')

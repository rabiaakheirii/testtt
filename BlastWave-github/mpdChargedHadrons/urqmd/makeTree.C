#include "stdio.h"
#include "fstream"
#include "iostream"
#include "TTree.h"
#include "TFile.h"

void makeTree( void )
{
    //string filename = "urqmd_output/full_evolution.txt";
    string filename = "urqmd_output/hyperCooperFray.txt";
    ifstream txtFile;
    txtFile.open(filename);

    float t, posX, posY, posZ;
    float e, P, T, QGPfraction;
    float fluidVelocityX, fluidVelocityY, fluidRapidity;
    float baryonDensity, muB, muS;
    float hypT, hypX, hypY, hypZ;
    
    TFile *f = new TFile("tree_hyper.root", "recreate"); //создаем файл
    TTree *tree = new TTree("tree", "Simple Tree"); //создаем дерево
    
    tree->Branch("t", &t, "t/F"); 
    tree->Branch("e", &e, "e/F"); 
    tree->Branch("T", &T, "T/F");
    tree->Branch("muB", &muB, "muB/F");

    int line = 0;
    while (true)
    {
        txtFile >> t >> posX >> posY >> posZ 
                >> e >> P >> T >> QGPfraction  
                >> fluidVelocityX >> fluidVelocityY >> fluidRapidity
                >> baryonDensity >> muB
                >> muS >> hypT >> hypX >> hypY >> hypZ;
                
        // cout << t << " " << posX << " " << posY << "  " << posZ << endl;;
        if (T > 0.01 && (hypX != 0 || hypY != 0 || hypZ != 0))
        {
            tree->Fill();
            // cout << "  T: " << T << " muB: " << muB << endl;
        }
        
        // if( line > 1000 ) break;
        if( txtFile.eof() ) break;
        line++;
    }
    
    txtFile.close();

    tree->Write();
    tree->Print();
}


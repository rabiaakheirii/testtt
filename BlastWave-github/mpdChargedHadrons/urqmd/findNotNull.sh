#awk '           
# f { if ($13 != "0.000000") print }
# /END_OF_HEADER/ { f=1 }
#' fort.21 > notnull.txt

awk '           
  f { if ($13 != "0.000000" && $13 != "-0.000000") print }
  /END_OF_HEADER/ { f=1 }
' fort.21 > hyperCooperFray.txt



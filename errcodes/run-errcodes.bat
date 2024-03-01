cd ..\rrLogit\src
..\..\errcodes\a.exe
R CMD BATCH icodes.R icodes.out
mv sysdata.rda ..\R\sysdata.rda
rm icodes.R
rm icodes.out
rm .RData
cd ..\..\errcodes


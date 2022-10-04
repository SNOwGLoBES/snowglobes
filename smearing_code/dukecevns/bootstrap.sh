#!/bin/sh
cd ..
if [ ! -d COHERENTProposal2018 ]; then
  git clone https://code.ornl.gov/COHERENT/COHERENTProposal2018.git
fi
cd COHERENTProposal2018
git pull
cd assumptions
rsync -rtphv eff gs jsonfiles qf *.cc *.root ../../dukecevns

cd ../../dukecevns
if [ ! -d json ]; then
  git clone https://github.com/nlohmann/json.git
else
  cd json; git pull; cd ..;
fi

make sns_rates

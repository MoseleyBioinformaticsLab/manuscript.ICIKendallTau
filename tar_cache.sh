cd /big_data/data/rmflight_ici-kendallt/
tar cvf manuscript.ICIKendallTau.targets_cache.tar _targets
tarsplitter -m split -i manuscript.ICIKendallTau.targets_cache.tar -o manuscript.ICIKendallTau.targets_cache.parts -p 3

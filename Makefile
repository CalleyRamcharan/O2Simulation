all: sim ana

sim: muon_gun_10GeV/generated_files/trddigits.root

ana: muon_gun_10GeV/generated_files/o2sim_HitsTRD.root
	cd $(dir $<); root -b -q -l '../../CheckHits.C()'

.PRECIOUS: muon_gun_10GeV/generated_files/o2sim_HitsTRD.root

muon_gun_10GeV/generated_files/o2sim_HitsTRD.root:
	mkdir -p $(dir $@)
	cd $(dir $@); o2-sim -m PIPE MAG TRD -n 10000 -g boxgen --configKeyValues 'BoxGun.pdg=13;BoxGun.eta[0]=-0.84;BoxGun.eta[1]=0.84;BoxGun.prange[0]=10;BoxGun.prange[1]=10'

muon_gun_10GeV/generated_files/trddigits.root: muon_gun_10GeV/generated_files/o2sim_HitsTRD.root
	cd $(dir $@); o2-sim-digitizer-workflow --onlyDet TRD

clean:
	rm -f muon_gun_10GeV/generated_files/*
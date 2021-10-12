NEVENTS ?= 1000

all: simp anap sime anae

simp: p_gun_2GeV/generated_files/trddigits.root

anap: p_gun_2GeV/generated_files/o2sim_HitsTRD.root
	cd $(dir $<); root -b -q -l '../../CheckHits.C()'

sime: e_gun_2GeV/generated_files/trddigits.root

anae: e_gun_2GeV/generated_files/o2sim_HitsTRD.root
	cd $(dir $<); root -b -q -l '../../CheckHits.C()'

.PRECIOUS: p_gun_2GeV/generated_files/o2sim_HitsTRD.root e_gun_2GeV/generated_files/o2sim_HitsTRD.root

p_gun_2GeV/generated_files/o2sim_HitsTRD.root:
	mkdir -p $(dir $@)
	cd $(dir $@); o2-sim -m PIPE MAG TRD -n 1000 -g boxgen --configKeyValues 'BoxGun.pdg=211;BoxGun.eta[0]=-0.84;BoxGun.eta[1]=0.84;BoxGun.prange[0]=2;BoxGun.prange[1]=2'

e_gun_2GeV/generated_files/o2sim_HitsTRD.root:
	mkdir -p $(dir $@)
	cd $(dir $@); o2-sim -m PIPE MAG TRD -n 1000 -g boxgen --configKeyValues 'BoxGun.pdg=11;BoxGun.eta[0]=-0.84;BoxGun.eta[1]=0.84;BoxGun.prange[0]=2;BoxGun.prange[1]=2'

p_gun_2GeV/generated_files/trddigits.root: p_gun_2GeV/generated_files/o2sim_HitsTRD.root
	cd $(dir $@); o2-sim-digitizer-workflow --onlyDet TRD

e_gun_2GeV/generated_files/trddigits.root: e_gun_2GeV/generated_files/o2sim_HitsTRD.root
	cd $(dir $@); o2-sim-digitizer-workflow --onlyDet TRD

clean:
	rm -f p_gun_2GeV/generated_files/*
	rm -f e_gun_2GeV/generated_files/*
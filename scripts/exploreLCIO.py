from pyLCIO import EVENT, UTIL, IOIMPL
import numpy as np
import matplotlib.pyplot as plt

inFile_noShield = "/scratch/jdervan/muon-collider-studies/output_sim/offsettest/neutron_evt10k_th90_ph0_1GeV_default.slcio" #change path to your desired input file

reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(inFile_noShield)

hitTimesAll = []
hitEnergiesAll = []

for evtNum, event in enumerate(reader):
    collectionEB    = event.getCollection("ECalBarrelCollection")
    encodingEB      = collectionEB.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
    decoderEB       = UTIL.BitField64(encodingEB)

    for i_hit, hit in enumerate(collectionEB):
        decoderEB.setValue((hit.getCellID0() & 0xFFFFFFFF) | (hit.getCellID1() << 32))
        layer = decoderEB['layer'].value()

        if layer > 0: #skip shielding layer
            if hit.getNMCContributions() > 0:
                hitTimes        = [hit.getTimeCont(i) for i in range(hit.getNMCContributions())]
                hitEnergies     = [hit.getEnergyCont(i) for i in range(hit.getNMCContributions())]

                if hitTimes: #check if hitTimes is not empty
                    hitTimesAll.extend(hitTimes)
                    hitEnergiesAll.extend(hitEnergies)

reader.close()


# -------- Plotting results -------- #
# Plot hit times, weighted by energy

fig, ax = plt.subplots(figsize=(8,6), dpi=100)
# Compute bins
# max_value = max(hitTimesAll)
max_value = 200
bin_width = max_value/50

bin_edges = np.arange(0, max_value + bin_width, bin_width)

ax.hist(hitTimesAll, weights=hitEnergiesAll, bins=bin_edges)
ax.set_yscale("log")
fig.savefig("/scratch/jdervan/muon-collider-studies/scripts/plots/offsettest/neutron_evt10k_th90_ph0_1GeV_default.png")

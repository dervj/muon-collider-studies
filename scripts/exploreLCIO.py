from pyLCIO import EVENT, UTIL, IOIMPL
import numpy as np
import matplotlib.pyplot as plt
import math
# from scipy.optimize import curve_fit
# from scipy.stats import chi2, moyal
# import scipy.special as sp




settings = {
    "filenames" : [
        "/scratch/jdervan/mucolltest/muon-collider-studies/output_sim/sim_neutrontest_500MeV_evt30k_ortho_bch0p05.slcio",
        # "/scratch/jdervan/mucolltest/muon-collider-studies/output_sim/sim_neutrontest_500MeV_evt30k_ortho_bch0p10.slcio",
        # "/scratch/jdervan/mucolltest/muon-collider-studies/output_sim/sim_neutrontest_500MeV_evt30k_ortho_bch0p25.slcio"
        # add additional comma-separated filenames here
    ],
    "labels"    : [
        "0.05*NIL BCH2",
        # "0.10*NIL BCH2",
        # "0.25*NIL BCH2",  


        # add legend labels for each additional input file here
    ]
}

hitTimesAll = [[] for _ in settings["filenames"]]
hitEnergiesAll = [[] for _ in settings["filenames"]]

for i, inFile in enumerate(settings["filenames"]):
    print("Processing",inFile,"...")
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(inFile)

    for evtNum, event in enumerate(reader):
        collectionEB    = event.getCollection("ECalBarrelCollection")
        encodingEB      = collectionEB.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoderEB       = UTIL.BitField64(encodingEB)

        for i_hit, hit in enumerate(collectionEB):
            decoderEB.setValue((hit.getCellID0() & 0xFFFFFFFF) | (hit.getCellID1() << 32))
            layer = decoderEB['layer'].value()

            if layer > 0: #skip shielding layer
                if hit.getNMCContributions() > 0:
                    hitTimes        = [hit.getTimeCont(j) for j in range(hit.getNMCContributions())]
                    hitEnergies     = [hit.getEnergyCont(j) for j in range(hit.getNMCContributions())]

                    if hitTimes: #check if hitTimes is not empty
                        hitTimesAll[i].extend(hitTimes)
                        hitEnergiesAll[i].extend(hitEnergies)

                        # if hasattr(hit, "getPosition"):
                        #     pos = hit.getPosition()
                        #     hitPositionsAll.append(hit.getPosition())
                        #     distance = math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
                        #     hitDistancesAll.append(distance)
                        # else:
                        #     pass

reader.close()

integrated_energies = []
for i in range(len(settings["filenames"])):
    energy_in_window = sum(
        energy for time, energy in zip(hitTimesAll[i], hitEnergiesAll[i])
        if 10 <= time <= 20
    )
    integrated_energies.append(energy_in_window)
    print(f"Integrated energy for file {settings['labels'][i]} (10 ns to 20 ns): {energy_in_window}")

fig, ax = plt.subplots(figsize=(8,6), dpi=100)

num_bins = 40
time_range = (10, 20)
bins = np.linspace(time_range[0], time_range[1], num_bins + 1)

for i in range(len(settings['filenames'])):
    ax.hist(
        hitTimesAll[i],
        bins=bins,
        range=time_range,
        weights=hitEnergiesAll[i],
        histtype="step",
        label=f"{settings['labels'][i]} (E_int = {integrated_energies[i]:.2f})"
    )
ax.set_xlabel("Hit time [ns]")
ax.set_ylabel("Energy-weighted counts (GeV)")
ax.set_title("500MeV neutron hit time distribution")
ax.legend()
ax.set_xlim(time_range)

fig.savefig("/scratch/jdervan/mucolltest/muon-collider-studies/scripts/500MeV_neutron_BCH2_scan.png")

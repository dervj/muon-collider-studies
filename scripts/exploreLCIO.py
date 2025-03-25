from pyLCIO import EVENT, UTIL, IOIMPL
import numpy as np
import matplotlib.pyplot as plt
import math

settings = {
    "filenames" : [
        "/scratch/jdervan/mucolltest/muon-collider-studies/output_sim/sim_neutrontest_500MeV_evt30k_ortho_bch0p05.slcio",
        "/scratch/jdervan/mucolltest/muon-collider-studies/output_sim/sim_neutrontest_500MeV_evt30k_ortho_bch0p10.slcio",
        "/scratch/jdervan/mucolltest/muon-collider-studies/output_sim/sim_neutrontest_500MeV_evt30k_ortho_bch0p25.slcio"
        # add additional comma-separated filenames here
    ],
    "labels"    : [
        "0.05*NIL BCH2",
        "0.10*NIL BCH2",
        "0.25*NIL BCH2",  
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

reader.close()

# –----- Calculating integrated energy in time window ------ 
integrated_energies         = []
integrated_energy_errors    = []
for i in range(len(settings["filenames"])):
    energies_window = [
        energy for time, energy in zip(hitTimesAll[i], hitEnergiesAll[i])
        if 10 <= time <= 20
    ]
    energy_in_window    = sum(energies_window)
    error_in_window     = math.sqrt(sum(energy**2 for energy in energies_window))
    integrated_energies.append(energy_in_window)
    integrated_energy_errors.append(error_in_window)
    print(f"Integrated energy for case {settings['labels'][i]} (10 ns to 20 ns): {energy_in_window:.2f} ± {error_in_window:.2f}")


# ----- Plotting the histograms with error bars -----
fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

num_bins    = 40
time_range  = (10, 20)
bins        = np.linspace(time_range[0], time_range[1], num_bins + 1)
bin_centers = (bins[:-1] + bins[1:]) / 2

for i in range(len(settings["filenames"])):
    hist, _         = np.histogram(hitTimesAll[i], bins=bins, range=time_range, weights=hitEnergiesAll[i])
    hist_err_sq, _  = np.histogram(hitTimesAll[i], bins=bins, range=time_range, weights=np.square(hitEnergiesAll[i]))
    hist_err        = np.sqrt(hist_err_sq)
    
    ax.hist(
        hitTimesAll[i],
        bins=bins,
        range=time_range,
        weights=hitEnergiesAll[i],
        histtype='step',
        label=f"{settings['labels'][i]} (E_int={integrated_energies[i]:.2f}±{integrated_energy_errors[i]:.2f})"
    )

    ax.errorbar(bin_centers, hist, yerr=hist_err, fmt='o', capsize=1, markersize='5')

ax.set_xlabel("Hit time [ns]")
ax.set_ylabel("Energy-weighted counts")
ax.set_title("500 MeV neutron simhit time distribution")
ax.legend()
ax.set_xlim(time_range)

fig.savefig("/scratch/jdervan/mucolltest/muon-collider-studies/scripts/500MeV_neutron_BCH2_scan.png")

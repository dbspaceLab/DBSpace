""" Higher order measures here"""
# Make a coherence generation function
def gen_coher(inpX, Fs=422, nfft=2**10, polyord=0):
    print("Starting a coherence run...")
    outPLV = nestdict()
    outCSD = nestdict()

    fvect = np.linspace(0, Fs / 2, nfft / 2 + 1)

    for chann_i in inpX.keys():
        print(chann_i)
        for chann_j in inpX.keys():
            csd_ensemble = np.zeros(
                (inpX[chann_i].shape[1], len(feat_order)), dtype=complex
            )
            plv = np.zeros((inpX[chann_i].shape[1], len(feat_order)))

            for seg in range(inpX[chann_i].shape[1]):
                # First we get the cross spectral density
                csd_out = sig.csd(
                    inpX[chann_i][:, seg], inpX[chann_j][:, seg], fs=Fs, nperseg=512
                )[1]

                # normalize the entire CSD for the total power in input signals
                norm_ms_csd = np.abs(csd_out) / np.sqrt(
                    l2_pow(inpX[chann_i][:, seg]) * l2_pow(inpX[chann_j][:, seg])
                )

                # Are we focusing on a band or doing the entire CSD?

                for bb, band in enumerate(feat_order):
                    # what are our bounds?
                    band_bounds = feat_dict[band]["param"]
                    band_idxs = np.where(
                        np.logical_and(fvect >= band_bounds[0], fvect <= band_bounds[1])
                    )
                    csd_ensemble[seg, bb] = np.median(csd_out[band_idxs])
                    plv[seg, bb] = np.max(norm_ms_csd[band_idxs])

                # Below brings in the entire csd, but this is dumb
                # csd_ensemble[seg] = csd_out

                # Compute the PLV

            # Here we find the median across segments
            outCSD[chann_i][chann_j] = csd_ensemble

            # Compute the normalized coherence/PLV
            outPLV[chann_i][chann_j] = plv
            # outPLV[chann_i][chann_j] = np.median(plv,axis=0)

            ## PLV abs EEG ->
            ## Coherence value

    return outCSD, outPLV

#%%
# BlackRock Methods
def load_or_file(fname, **kwargs):
    # nsx_file = NsxFile(fname)

    arg_list = ["elec_ids", "start_time_s", "data_time_s", "downsample", "plot_chan"]

    for aa in arg_list:
        print(arg_list[aa])

    # Extract data - note: data will be returned based on *SORTED* elec_ids, see cont_data['elec_ids']
    cont_data = nsx_file.getdata(elec_ids, start_time_s, data_time_s, downsample)

    # Close the nsx file now that all data is out
    nsx_file.close()

    return cont_data

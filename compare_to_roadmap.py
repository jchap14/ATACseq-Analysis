def compare_to_roadmap(bw_file, regions_file, reg2map_file,
                       metadata, output_prefix):
    '''
    Takes a bigwig file and signal file, gets the bwAverageOverBed,
    then compares that signal with the signal in the Roadmap
    regions
    '''

    out_file = '{0}.signal'.format(output_prefix)

    # First get the signal vals for the peak regions
    # remember to use a UCSC formatted bed file for regions
    bw_average_over_bed = 'bigWigAverageOverBed {0} {1} {2}'.format(
                            bw_file, regions_file, out_file)
    logging.info(bw_average_over_bed)
    os.system(bw_average_over_bed)

    # Read the file back in
    sample_data = pd.read_table(out_file, header=None)
    sample_mean0_col = np.array(sample_data.iloc[:, 5])

    # Then, calculate correlations with all other Roadmap samples and rank
    # the correlations
    roadmap_signals = pd.read_table(reg2map_file, compression='gzip')
    (nrow, ncol) = roadmap_signals.shape

    results = pd.DataFrame(columns=('eid', 'corr'))
    for i in range(ncol):
        # Slice, run correlation
        roadmap_i = roadmap_signals.iloc[:, i]
        spearman_corr = scipy.stats.spearmanr(np.array(roadmap_i),
                                              sample_mean0_col)
        results.loc[i] = [roadmap_i.name, spearman_corr[0]]
        logging.info('{0}\t{1}'.format(roadmap_i.name, spearman_corr))

    # Read in metadata to make the chart more understandable
    metadata = pd.read_table(metadata)
    metadata.columns = ['eid', 'mnemonic']

    merged = pd.merge(metadata, results, on='eid')

    sorted_results = merged.sort('corr', ascending=True)

    # Plot results
    pos = np.array(range(ncol)) + 0.5
    fig = plt.figure()
    plt.barh(pos, sorted_results['corr'], align='center', height=1.0)
    plt.yticks(pos, sorted_results['mnemonic'].tolist(), fontsize=7)
    plt.xlabel('Spearmans correlation')
    plt.title('Signal correlation to Roadmap DNase')
    plt.axis('tight')
    ax = plt.axes()
    ax.yaxis.set_ticks_position('none')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png', bbox_inches='tight')
    fig.savefig('test.png', format='png', bbox_inches='tight')

    return b64encode(plot_img.getvalue())
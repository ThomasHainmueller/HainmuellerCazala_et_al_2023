import sima

def correct_sbx(fname,max_displacement=[30,30], num_states_retained=50,
                correction_channel='red', bigfile=False, mode = 1, nlines = 8):
    """
    Wrapper function to use sima motion correction from Matlab
    """

    # Generate the 2D-HMM by which all files will be corrected
    if mode == 1:
    	markov2D = sima.motion.HiddenMarkov2D(
            granularity='row',num_states_retained=num_states_retained,
            max_displacement=max_displacement, n_processes=1);
    
    # Alternative 2D-HMM for IN data with sparse labelling
    elif mode == 2:
        markov2D = sima.motion.HiddenMarkov2D(
            granularity=('row',nlines),num_states_retained=num_states_retained,
            max_displacement=max_displacement, n_processes=1);
    
    sequences = [sima.Sequence.create('SBX',fname)]

    dataset = markov2D.correct(
        sequences, fname.replace('.sbx','.sima'), channel_names=['green','red'],
        correction_channels=[correction_channel])

    if bigfile:
        nframes = len(dataset.sequences[0])
        return iter(dataset.sequences[0]), nframes
    else:
        return dataset.sequences[0].export_matlab()

def load_simafile(simaname, bigfile=False):
    """
    Wrapper function to load a motion-corrected sima-sequence from matlab.
    """
    dataset = sima.ImagingDataset.load(simaname)

    if bigfile:
        nframes = len(dataset.sequences[0])
        return iter(dataset.sequences[0]), nframes
    else:
        return dataset.sequences[0].export_matlab()

def create_iterator(simaname):
    """
    Create iterator for large sima files that can't be exported at once
    """
    dataset = sima.ImagingDataset.load(simaname)
    nframes = len(dataset.sequences[0])
    
    return iter(dataset.sequences[0]), nframes

from collections import defaultdict
from typing import List

from . import ConsensusMap, ConsensusFeature, FeatureMap, Feature, MSExperiment, PeakMap, PeptideIdentification, ControlledVocabulary, File, IonSource

import pandas as pd
import numpy as np

class ConsensusMapDF(ConsensusMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_intensity_df(self):
        """Generates a pandas DataFrame with feature intensities from each sample in long format (over files).

        For labelled analyses channel intensities will be in one row, therefore resulting in a semi-long/block format.
        Resulting DataFrame can be joined with result from get_metadata_df by their index 'id'.

        Returns:
        pandas.DataFrame: intensity DataFrame
        """
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()  # type: dict[int, ColumnHeader]

        labels = list(set([header.label for header in filemeta.values()]))
        files = list(set([header.filename for header in filemeta.values()]))
        label_to_idx = {k: v for v, k in enumerate(labels)}
        file_to_idx = {k: v for v, k in enumerate(files)}

        def gen(cmap: ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        if not labelfree:

            def extract_row_blocks_channel_wide_file_long(f: ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                filerows = defaultdict(lambda: [0] * len(labels))
                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row = filerows[header.filename]
                    row[label_to_idx[header.label]] = fh.getIntensity()
                return (f.getUniqueId(), filerows)

            def extract_rows_channel_wide_file_long(f: ConsensusFeature):
                uniqueid, rowdict = extract_row_blocks_channel_wide_file_long(f)
                for file, row in rowdict.items():
                    row.append(file)
                    yield tuple([uniqueid] + row)

            if len(labels) == 1:
                labels[0] = "intensity"

            dtypes = [('id', np.dtype('uint64'))] + list(zip(labels, ['f'] * len(labels)))
            dtypes.append(('file', 'U300'))

            intyarr = np.fromiter(iter=gen(self, extract_rows_channel_wide_file_long), dtype=dtypes, count=self.size())

            return pd.DataFrame(intyarr).set_index('id')

        else:
            # Specialized for LabelFree which has to have only one channel
            def extract_row_blocks_channel_long_file_wide_LF(f: ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                row = [0.] * len(files)

                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row[file_to_idx[header.filename]] = fh.getIntensity()

                yield tuple([f.getUniqueId()] + row)

            dtypes = [('id', np.dtype('uint64'))] + list(zip(files, ['f'] * len(files)))

            intyarr = np.fromiter(iter=gen(self, extract_row_blocks_channel_long_file_wide_LF), dtype=dtypes, count=self.size())

            return pd.DataFrame(intyarr).set_index('id')

    def get_metadata_df(self):
        """Generates a pandas DataFrame with feature meta data (sequence, charge, mz, RT, quality).

        Resulting DataFrame can be joined with result from get_intensity_df by their index 'id'.

        Returns:
        pandas.DataFrame: DataFrame with metadata for each feature (such as: best identified sequence, charge, centroid RT/mz, fitting quality)
        """

        def gen(cmap: ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        def extract_meta_data(f: ConsensusFeature):
            pep = f.getPeptideIdentifications()  # type: list[PeptideIdentification]

            if len(pep) != 0:
                hits = pep[0].getHits()

                if len(hits) != 0:
                    besthit = hits[0]  # type: PeptideHit
                    yield f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()
                
                else:
                    yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()
            
            else:
                yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()

        cnt = self.size()

        mddtypes = [('id', np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'),
                    ('RT', np.dtype('double')), ('mz', np.dtype('double')), ('quality', 'f')]

        mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return pd.DataFrame(mdarr).set_index('id')
    
    def get_df(self):
        """Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        Returns:
        pandas.DataFrame: meta data and intensity DataFrame
        """
        return pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)

ConsensusMap = ConsensusMapDF

# TODO tell the advanced user that they could change this, in case they have different needs.
# TODO check if type could be inferred in the first pass
# TODO check if max. string lengths could be inferred in the first pass and how this affects runtime
# TODO check how runtime is affected if we use np.append instead of np.fromiter and use np.dyte = object for strings
common_meta_value_types = {
    b'label': 'U50',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'U100',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'f', # TODO this is actually a DoubleList. Think about what to do here. For np.fromiter we would need to set the length of the array.
    b'Group': 'U50',
    b'is_ungrouped_monoisotopic': 'i' # TODO this sounds very boolean to me
}

class FeatureMapDF(FeatureMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __get_prot_id_filename_from_pep_id(self, pep_id):
        """Gets the primary MS run path of the ProteinIdentification linked with the given PeptideIdentification.

        Parameters:
        pep_id: PeptideIdentification

        Returns:
        str: primary MS run path (filename) of the ProteinIdentification with the same identifier as the given PeptideIdentification
        """
        for prot in self.getProteinIdentifications():
            if prot.getIdentifier() == pep_id.getIdentifier():
                filenames = []
                prot.getPrimaryMSRunPath(filenames)
                if filenames and filenames[0] != '':
                    return filenames[0]
        return 'unknown'
    
    # meta_values = None (default), 'all' or list of meta value names
    def get_df(self, meta_values = None, export_peptide_identifications = True):
        """Generates a pandas DataFrame with information contained in the FeatureMap.

        Optionally the feature meta values and information for the assigned PeptideHit can be exported.

        Parameters:
        meta_values: meta values to include (None, [custom list of meta value names] or 'all')

        export_peptide_identifications (bool): export sequence and score for best PeptideHit assigned to a feature.
        Additionally the ID_filename (file name of the corresponding ProteinIdentification) and the ID_native_id 
        (spectrum ID of the corresponding Feature) are exported. They are also annotated as meta values when 
        collecting all assigned PeptideIdentifications from a FeatureMap with FeatureMap.get_assigned_peptide_identifications().
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """
        # get all possible meta value keys in a set
        if meta_values == 'all':
            meta_values = set()
            for f in self:
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values.add(m)

        elif not meta_values: # if None, set to empty list
            meta_values = []
        
        def gen(fmap: FeatureMap, fun):
            for f in fmap:
                yield from fun(f)

        def extract_meta_data(f: Feature):
            """Extracts feature meta data.
            
            Extracts information from a given feature with the requested meta values and, if requested,
            the sequence, score and ID_filename (primary MS run path of the linked ProteinIdentification)
            of the best PeptideHit (first) assigned to that feature.

            Parameters:
            f (Feature): feature from which to extract the meta data

            Yields:
            tuple: tuple containing feature information, peptide information (optional) and meta values (optional)
            """
            pep = f.getPeptideIdentifications()  # type: list[PeptideIdentification]
            bb = f.getConvexHull().getBoundingBox2D()
                
            vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.nan for m in meta_values]
            
            if export_peptide_identifications:
                if len(pep) > 0:
                    ID_filename = self.__get_prot_id_filename_from_pep_id(pep[0])
                    hits = pep[0].getHits()
                    if len(hits) > 0:
                        besthit = hits[0]
                        pep_values = (besthit.getSequence().toString(), besthit.getScore(), ID_filename, f.getMetaValue('spectrum_native_id'))
                    else:
                        pep_values = (None, None, ID_filename, f.getMetaValue('spectrum_native_id'))
                else:
                    pep_values = (None, None, None, None)
            else:
                pep_values = ()

            yield tuple([f.getUniqueId()]) + pep_values + (f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals)

        cnt = self.size()

        mddtypes = [('feature_id', 'U100')]
        if export_peptide_identifications:
            mddtypes += [('peptide_sequence', 'U200'), ('peptide_score', 'f'), ('ID_filename', 'U100'), ('ID_native_id', 'U100')]
        mddtypes += [('charge', 'i4'), ('RT', np.dtype('double')), ('mz', np.dtype('double')), ('RTstart', np.dtype('double')), ('RTend', np.dtype('double')),
                    ('MZstart', np.dtype('double')), ('MZend', np.dtype('double')), ('quality', 'f'), ('intensity', 'f')]
        
        for meta_value in meta_values:
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
            else:
                mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return pd.DataFrame(mdarr).set_index('feature_id')

    def get_assigned_peptide_identifications(self):
        """Generates a list with peptide identifications assigned to a feature.

        Adds 'ID_native_id' (feature spectrum id), 'ID_filename' (primary MS run path of corresponding ProteinIdentification)
        and 'feature_id' (unique ID of corresponding Feature) as meta values to the peptide hits.
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])

        Returns:
        [PeptideIdentification]: list of PeptideIdentification objects
        """
        result = []
        for f in self:
            for pep in f.getPeptideIdentifications():
                hits = []
                for hit in pep.getHits():
                    hit.setMetaValue('feature_id', str(f.getUniqueId()))
                    hit.setMetaValue('ID_filename', self.__get_prot_id_filename_from_pep_id(pep))
                    if f.metaValueExists('spectrum_native_id'):
                        hit.setMetaValue('ID_native_id', f.getMetaValue('spectrum_native_id'))
                    else:
                        hit.setMetaValue('ID_native_id', 'unknown')
                    hits.append(hit)
                pep.setHits(hits)
                result.append(pep)
        return result

FeatureMap = FeatureMapDF


class MSExperimentDF(MSExperiment):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, long : bool = False):
        """Generates a pandas DataFrame with all peaks in the MSExperiment

        Parameters:
        long: set to True if you want to have a long/expanded/melted dataframe with one row per peak. Faster but
            replicated RT information. If False, returns rows in the style: rt, np.array(mz), np.array(int)
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """
        if long:
            cols = ["RT", "mz", "inty"]
            self.updateRanges()
            spectraarrs2d = self.get2DPeakDataLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ())
            return pd.DataFrame(dict(zip(cols, spectraarrs2d)))

        cols = ["RT", "mzarray", "intarray"]

        return pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in self), columns=cols)

    def get_massql_df(self):
        """Exports data from MSExperiment to pandas DataFrames to be used with MassQL.

        The Python module massql allows queries in mass spectrometry data (MS1 and MS2
        data frames) in a SQL like fashion (https://github.com/mwang87/MassQueryLanguage).
        
        Both dataframes contain the columns:
        'i': intensity of a peak
        'i_norm': intensity normalized by the maximun intensity in the spectrum
        'i_tic_norm': intensity normalized by the sum of intensities (TIC) in the spectrum
        'mz': mass to charge of a peak
        'scan': number of the spectrum
        'rt': retention time of the spectrum
        'polarity': ion mode of the spectrum as integer value (positive: 1, negative: 2)
        
        The MS2 dataframe contains additional columns:
        'precmz': mass to charge of the precursor ion
        'ms1scan': number of the corresponding MS1 spectrum
        'charge': charge of the precursor ion

        Returns:
        ms1_df (pandas.DataFrame): peak data of MS1 spectra
        ms2_df (pandas.DataFrame): peak data of MS2 spectra with precursor information
        """
        self.updateRanges()

        def _get_polarity(spec):
            '''Returns polarity as an integer value for the massql dataframe.
            
            According to massql positive polarity is represented by 1 and negative by 2.

            Parameters:
            spec (MSSpectrum): the spectrum to extract polarity

            Returns:
            int: polarity as int value according to massql specification
            '''
            polarity = spec.getInstrumentSettings().getPolarity()
            if polarity == IonSource.Polarity.POLNULL:
                return 0
            elif polarity == IonSource.Polarity.POSITIVE:
                return 1
            elif polarity == IonSource.Polarity.NEGATIVE:
                return 2

        def _get_spec_arrays(mslevel):
            '''Get spectrum data as a matrix.

            Generator yields peak data from each spectrum (with specified MS level) as a numpy.ndarray.
            Normalized intensity values are calculated and the placeholder values replaced. For 'i_norm' and
            'i_tic_norm' the intensity values are divided by the maximum intensity value in the spectrum and 
            the sum of intensity values, respectively.

            Parameters:
            mslevel (int): only spectra with the given MS level will be considered

            Yields:
            np.ndarray: 2D array with peak data (rows) from each spectrum
            '''
            for scan_num, spec in enumerate(self):
                if spec.getMSLevel() == mslevel:
                    mz, inty = spec.get_peaks()
                    # data for both DataFrames: i, i_norm, i_tic_norm, mz, scan, rt, polarity
                    data = (inty, inty/np.amax(inty, initial=0), inty/np.sum(inty), mz, scan_num + 1, spec.getRT()/60, _get_polarity(spec))
                    cols = 7
                    if mslevel == 2:
                        cols = 10
                        # data for MS2 only: precmz, ms1scan, charge
                        # set fallback values if no precursor is annotated (-1)
                        if spec.getPrecursors():
                            data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num)+1, spec.getPrecursors()[0].getCharge())
                        else:
                            data += (-1, -1, -1)
                    # create empty ndarr with shape according to MS level
                    ndarr = np.empty(shape=(spec.size(), cols))
                    # set column values
                    for i in range(cols):
                        ndarr[:,i] = data[i]
                    yield ndarr

        # create DataFrame for MS1 and MS2 with according column names and data types
        # if there are no spectra of given MS level return an empty DataFrame
        dtypes = {'i': 'float32', 'i_norm': 'float32', 'i_tic_norm': 'float32', 'mz': 'float64', 'scan': 'int32', 'rt': 'float32', 'polarity': 'int32'}
        if 1 in self.getMSLevels():
            ms1_df = pd.DataFrame(np.concatenate(list(_get_spec_arrays(1)), axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms1_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        dtypes = dict(dtypes, **{'precmz': 'float64', 'ms1scan': 'int32', 'charge': 'int32'})
        if 2 in self.getMSLevels():
            ms2_df = pd.DataFrame(np.concatenate(list(_get_spec_arrays(2)), axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms2_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        return ms1_df, ms2_df

MSExperiment = MSExperimentDF
PeakMap = MSExperimentDF


# TODO think about the best way for such top-level function. IMHO in python, encapsulation in a stateless class in unnecessary.
#   We should probably not just import this whole submodule without prefix.
def peptide_identifications_to_df(peps: List[PeptideIdentification], decode_ontology : bool = True,
                                  default_missing_values: dict = {bool: False, int: -9999, float: np.nan, str: ''},
                                  export_unidentified : bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.
    Parameters:
    peps (List[PeptideIdentification]): list of PeptideIdentification objects
    decode_ontology (bool): decode meta value names
    default_missing_values: default value for missing values for each data type
    export_unidentified: export PeptideIdentifications without PeptideHit
    Returns:
    pandas.DataFrame: peptide identifications in a DataFrame
    """
    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}

    # filter out PeptideIdentifications without PeptideHits if export_unidentified == False
    count = len(peps)
    if not export_unidentified:
        count = sum(len(pep.getHits()) > 0 for pep in peps)

    # get all possible metavalues
    metavals = []
    types = []
    mainscorename = "score"
    for pep in peps:
        hits = pep.getHits()
        if not len(hits) == 0:
            mvs = []
            hits[0].getKeys(mvs)
            metavals += mvs
            mainscorename = pep.getScoreType()

    metavals = list(set(metavals))

    # get type of all metavalues
    for k in metavals:
        if k == b"target_decoy":
            types.append('?')
        else:
            for p in peps:
                hits = p.getHits()
                if not len(hits) == 0:
                    mv = hits[0].getMetaValue(k)
                    types.append(switchDict[type(mv)])
                    break

    # get default value for each type in types to append if there are no hits in a PeptideIdentification
    def get_key(val):
        for key, value in switchDict.items():
            if val == value:
                return key
    dmv = [default_missing_values[get_key(t)] for t in types]

    decodedMVs = [m.decode("utf-8") for m in metavals]
    if decode_ontology:
        cv = ControlledVocabulary()
        cv.loadFromOBO("psims", File.getOpenMSDataPath() + "/CV/psi-ms.obo")
        clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
    else:
        clearMVs = decodedMVs
        
    clearcols = ["id", "RT", "mz", mainscorename, "charge", "protein_accession", "start", "end"] + clearMVs
    coltypes = ['U100', 'f', 'f', 'f', 'i','U1000', 'U1000', 'U1000'] + types
    dt = list(zip(clearcols, coltypes))

    def extract(pep):
        hits = pep.getHits()
        if not hits:
            if export_unidentified:
                return (pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), default_missing_values[float], default_missing_values[int],
                        default_missing_values[str], default_missing_values[str], default_missing_values[str], *dmv)
            else:
                return

        besthit = hits[0]
        ret = [pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), besthit.getScore(), besthit.getCharge()] 
        # add accession, start and end positions of peptide evidences as comma separated str (like in mzTab)
        evs = besthit.getPeptideEvidences()
        ret += [','.join(v) if v else default_missing_values[str] for v in ([e.getProteinAccession() for e in evs],
                                                                            [str(e.getStart()) for e in evs],
                                                                            [str(e.getEnd()) for e in evs])]
        for k in metavals:
            if besthit.metaValueExists(k):
                val = besthit.getMetaValue(k)
                if k == b"target_decoy":
                    if val[0] == 't':
                        ret.append(True)
                    else:
                        ret.append(False)
                else:
                    ret.append(val)
            else:
                ret.append(default_missing_values[type(val)])
        return tuple(ret)

    return pd.DataFrame(np.fromiter((extract(pep) for pep in peps), dtype=dt, count=count))

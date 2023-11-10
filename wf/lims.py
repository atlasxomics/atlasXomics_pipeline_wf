import csv

from slims.criteria import equals
from slims.slims import Slims

from wf.creds import username, password


NA = ['na', '#N/A', '=NA()', 'NA', 'nan', None]


def csv_to_dict(path):
    with open(path) as f:
        reader = csv.reader(f)
        return dict(zip(reader.__next__(), reader.__next__()))


def get_pk(cntn_id, slims):
    return slims.fetch('Content', equals('cntn_id', cntn_id))[0].pk()


def push_result(payload, slims):
    return slims.add("Result", payload)


def slims_init(username=username, password=password):
    return Slims(
        "slims", "https://slims.atlasxomics.com/slimsrest", username, password
    )


mapping = {
    'Genome': 'rslt_cf_refGenome',
    'Pipeline version': 'rslt_cf_pipelineVersion',
    'Fraction aligned reads ': 'rslt_cf_confidentlyMappedReadPairs',
    'FRIP': 'rslt_cf_fractionOfHighQualityFragmentsOrlapPe',
    'Number of peaks': 'rslt_cf_numberOfPeaks',
    'Fraction duplicate reads': 'rslt_cf_percentDuplicates',
    'Chromap input read pairs': 'rslt_cf_sequencedReadPairs1',
    'TSS_enrichment': 'rslt_cf_tssEnrichmentScore',
    'Fraction unaligned reads': 'rslt_cf_unmappedReadPairs',
    'Fraction reads with valid barcode': 'rslt_cf_validBarcodes'
}

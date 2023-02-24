from ..session import Component
from ..referencebuilder import ReferenceBuilder
from pathlib import Path
from hashlib import md5
import pytest
import shutil

@pytest.fixture
def builder_defec(tmpdir, datadir):
    return  ReferenceBuilder(Component(ref_dir=Path(tmpdir)), input_fasta = datadir/ "defect_catalogue.fasta.gz")


@pytest.fixture
def builder(tmpdir, datadir):
    return  ReferenceBuilder(Component(ref_dir=Path(tmpdir), ref_name="test"), input_fasta = datadir/ "test_catalogue.fasta.gz")

# def get_read_reference(builder_defec):
#     yield from builder_defec.read_reference()

# @pytest.fixture(params=get_read_reference())
# def line_fixture(request):
#     return request.param

# @pytest.mark.parametrize(
#     ("length", "sequence"),
#     (
#         (80, "ATGAAAATGAACCTGCAAAAGGACATGTTTGATCGCAAATTGCGATACAAGATGTACAAAGATGGTAAAAAGTGGGTGTT"),
#         (56, "TGCCAGCATGGCGACTTTGTCTTTGATTGGGGCTTTTTTGGGCGGTGGTAGCGCCC"),
#         (74, "ATGAAAATGAACCTGCAAAAGGACGATCGCAAATTGCGATACAAGATGTACAAAGATGGTAAAAAGTGGGTGTT"),
#     ),
# )
def test_read_reference(builder_defec):
    res =  (
        (80, "ATGAAAATGAACCTGCAAAAGGACATGTTTGATCGCAAATTGCGATACAAGATGTACAAAGATGGTAAAAAGTGGGTGTT"),
        (56, "TGCCAGCATGGCGACTTTGTCTTTGATTGGGGCTTTTTTGGGCGGTGGTAGCGCCC"),
        (74, "ATGAAAATGAACCTGCAAAAGGACGATCGCAAATTGCGATACAAGATGTACAAAGATGGTAAAAAGTGGGTGTT"),
    )
    for i, (length, sequence) in enumerate(builder_defec.read_reference()):
        assert length == res[i][0]
        assert sequence == res[i][1]

@pytest.mark.parametrize(
    ("annotation_md5", "fasta_md5"),
    (
        pytest.param("b8d7e70fe88abbf1cf1746e5a02439ea", "2912b682a8e7554025cc5feadd641570", id="Accurate output"),
    )
)
def test_create_reference(builder, annotation_md5, fasta_md5):
    builder.create_reference()
    with open(builder.output_annotation_file, "rb") as output_annotation_file:
        assert md5(output_annotation_file.read()).hexdigest() == annotation_md5
    with open(builder.output_fasta_file, "rb") as output_fasta_file:
        assert md5(output_fasta_file.read()).hexdigest() == fasta_md5

    # with pytest.raises(ValueError):
    #     Line that raise a value error

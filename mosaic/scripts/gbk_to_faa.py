from Bio import SeqIO
import sys
gbk_filename=sys.argv[1]
faa_filename=sys.argv[2]
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")
print(input_handle)
for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            print(seq_feature)
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()

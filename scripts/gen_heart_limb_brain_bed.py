from Bio import SeqIO

def description_to_bedline(description, outputFile):
    description = description.split("|")[1]
    chromName = description.split(":")[0]

    chromRange = description.split(":")[1]
    chromStart = chromRange.split("-")[0]
    chromEnd = chromRange.split("-")[1]

    bedLine = [chromName, chromStart, chromEnd]
    bedLine = '\t'.join(bedLine) + "\n"
    
    outputFile.write(bedLine)

heartFile = open("data/bed/heart.bed", 'w')
limbFile = open("data/bed/limb.bed", 'w')
brainFile = open("data/bed/brain.bed", 'w')
comboFile = open("data/bed/combo.bed", 'w')

FH_f = open("data/fasta/hm_annotated.fa")
human_fasta_seq = SeqIO.parse(FH_f, 'fasta')
seqs = []

heart = "heart"
limb = "limb"
brain = "brain"

for x in human_fasta_seq:
    description = x.description
    if heart in description and limb not in description and brain not in description:
        description_to_bedline(description, heartFile)
    elif limb in description and heart not in description and brain not in description:
        description_to_bedline(description, limbFile)
    elif brain in description and limb not in description and heart not in description:
        description_to_bedline(description, brainFile)
    elif brain in description or limb in description or heart in description:
        description_to_bedline(description, comboFile)


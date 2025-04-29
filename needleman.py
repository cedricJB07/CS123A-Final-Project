#read first 150 sequences from the FASTQ file
def get_sequences(filename):
    sequences = []
    with open(filename) as file:
        lines = file.readlines()
        #every 4 lines, line 2 (index 1) is the sequence based on FastQ structure
        for i in range(1, len(lines), 4):  
            sequences.append(lines[i].strip())
            if len(sequences) == 150:
                break
    print(f"Extracted {len(sequences)} sequences from {filename}")
    return sequences


#align two sequences using Needleman-Wunsch
def needleman_wunsch(seq1, seq2):

    #scoring
    match = 1
    mismatch = -1
    gap = -2

    #create the scoring matrix
    m, n = len(seq1), len(seq2)
    score = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    #initialize the first row and first column
    for i in range(m + 1):
        score[i][0] = i * gap
    for j in range(n + 1):
        score[0][j] = j * gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            #diagonal move (match or mismatch)
            if seq1[i - 1] == seq2[j - 1]:
                diag = score[i - 1][j - 1] + match
            else:
                diag = score[i - 1][j - 1] + mismatch
            #up move 
            up = score[i - 1][j] + gap
            #left move
            left = score[i][j - 1] + gap

            #pick the highest score among diagonal,up,left
            score[i][j] = max(diag, up, left)

            
    #return the bottom right most score after the matrix is filled in
    return score[m][n]



#files
ground_file = "GLDS-621_rna-seq_Bulk_072022_1G_R1_raw.fastq"
microG_file = "GLDS-621_rna-seq_Bulk_072022_uG_R1_raw.fastq"

ground_reads = get_sequences(ground_file)
space_reads = get_sequences(microG_file)


# Align all 150 reads and calculate the average alignment score
total_score = 0

print("Aligning all 150 Ground vs Microgravity reads: ")


for i in range(150):
    score = needleman_wunsch(ground_reads[i], space_reads[i])
    print(f"Alignment score for Read {i+1}: {score}")
    total_score += score

average_score = total_score / 150


print("\nFinal average alignment score across all 150 reads:", average_score)

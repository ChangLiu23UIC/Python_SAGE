def decoy_generation(fasta_file:str):
    """
    Generates decoy for the FDR control
    :param fasta_file:
    :return:
    """
    with open(fasta_file, "r") as file:
        with open(fasta_file.split(".")[0]+"_rev.fasta", "w") as writing:
            file_content = file.read().split("\n>")
            for entry in file_content:
                lines = entry.strip().split('\n')
                header = ">" + lines[0] + "\n"
                seq = "".join(lines[1:])
                seq_n = seq + "\n"
                rev_header = ">rev_" + lines[0] + "\n"
                rev_seq = seq[::-1] + "\n"
                writing.write(header)
                writing.write(seq_n)
                writing.write(rev_header)
                writing.write(rev_seq)
        writing.close()
        return file_content

if __name__ == '__main__':
    dd = decoy_generation("human_1.fasta")
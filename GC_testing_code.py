import os

directory = "332_genome_assemblies"

def count_characters(file_path):
    try:
        with open(file_path, 'r') as file:
            atcg_count = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
            lower_atcg_count = {'a':0, 't':0, 'c':0, 'g':0}
            char_count = 0
            for line in file:
                if line.startswith('>'):
                    continue
                else:
                    char_count += len(line)
                    atcg_count['A'] += line.count('A')
                    lower_atcg_count['a'] += line.count('a')

                    atcg_count['T'] += line.count('T')
                    lower_atcg_count['t'] += line.count('t')

                    atcg_count['C'] += line.count('C')
                    lower_atcg_count['c'] += line.count('c')

                    atcg_count['G'] += line.count('G')
                    lower_atcg_count['g'] += line.count('g')

                    atcg_count['N'] += line.count('N')
            atcg_total = sum(atcg_count.values()) + sum(lower_atcg_count.values()) - atcg_count['N']
            n_atcg_total = sum(atcg_count.values())
            lower_atcg_total = sum(lower_atcg_count.values())

            print(f"A: {atcg_count['A']} occurrences")
            print(f"T: {atcg_count['T']} occurrences")
            print(f"C: {atcg_count['C']} occurrences")
            print(f"G: {atcg_count['G']} occurrences\n")

            print("Total ATCG count:", atcg_total)
            print("Total ATCG count WITH N:", n_atcg_total,"\n")
            print("Total LOWERCASE ATCG count:", lower_atcg_total,"\n")
            gc_content = ((atcg_count['G']+lower_atcg_count['g'] + atcg_count['C']+lower_atcg_count['c']) / atcg_total) * 100
            print("GC content:", gc_content, "%\n")
            print(f"N: {atcg_count['N']} occurrences")

            return char_count, atcg_count, gc_content  # Return char_count, atcg_count, and gc_content
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")

def process_files(folder_path):
    with open("TEST.csv", 'a') as output_file:
        output_file.write("File,Total Characters,A,T,C,G,GC Content,N\n")

        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)

            if os.path.isfile(file_path) and file_name.endswith('.fas'):
                print("Processing file:", file_name)
                result, atcg_count, gc_content = count_characters(file_path)
                output_file.write(f"{file_name},{result},{atcg_count['A']},{atcg_count['T']},{atcg_count['C']},{atcg_count['G']},{gc_content},{atcg_count['N']}\n")

file_path = directory
process_files(file_path)

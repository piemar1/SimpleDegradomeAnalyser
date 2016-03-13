# ! /usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Marcin Pieczyński'

import datetime


    # plik = "/media/marcin/01D068D06AB0C3C0/SimpleDegradomeSequenseAnalyser/root_S10_L001_R1_001.fastq"   # ubuntu
    # plik_output = "/media/marcin/01D068D06AB0C3C0/SimpleDegradomeSequenseAnalyser/NEWOUTPUT1.txt"   # ubuntu

    # plik = "D:\SimpleDegradomeSequenseAnalyser\\root_S10_L001_R1_001.fastq"                              # win
    # plik_output = "D:\SimpleDegradomeSequenseAnalyser\\NEWOUTPUT2.txt"                              # win

    # adapter_3 = "ATGGG"
    # adapter_5 = "TGGAATTCTC"
    # sorting_way = "L"


class DegradomeSequenceAnalyser(object,):
    def __init__(self):
        # self.adapter_3 = "ATGGG"
        # self.adapter_5 = "TGGAATTCTC"
        self.number_of_parts = 16           # number of data fragments do not change !!

        self.clustered_data = []
        self.seq_3_5_cuted = []

        self.clustered_collected_A = []
        self.clustered_collected_B = []
        self.clustered_collected_C = []
        self.clustered_collected_D = []

        self.input_file_readlines = None

    def open_input_file(self, input_file_path):
        self.input_file_readlines = open(input_file_path).readlines()


    def filtr_cut_3_adaptors(self, adaptor3):

        def filtr_3_adaptors(elem, adaptor3):
            if adaptor3 in elem:
                return True
            return False

        def cut_3_adaptors(elem):
            cut = elem.index(adaptor3)
            return elem[cut:].strip()

        self.seq_3_5_cuted = map(cut_3_adaptors, [elem for elem in self.input_file_readlines
                                                    if filtr_3_adaptors(elem, adaptor3)])
        self.input_file_readlines = [elem for elem in self.seq_3_5_cuted if elem]


    def filtr_cut_5_adaptors(self, adaptor5):

        def filtr_5_adaptors(elem, adaptor5):
            if adaptor5 in elem:
                return True
            return False

        def cut_5_adaptors(elem):
            cut = elem.index(adaptor5)
            return elem[:cut].strip()

        self.seq_3_5_cuted = map(cut_5_adaptors, [elem for elem in self.input_file_readlines
                                                    if filtr_5_adaptors(elem, adaptor5)])

        self.input_file_readlines = [elem for elem in self.seq_3_5_cuted if elem]

    # def filtr_cut_3_5_adaptors(self, adaptor3, adaptor5):
    #
    #     def filtr_3_5_adaptors(elem, adaptor3, adaptor5):
    #         if adaptor5 in elem and elem[0:5] == adaptor3:
    #             return True
    #         return False
    #     def cut_3_5_adaptors(elem):
    #         cut = elem.index(adaptor5)
    #         return elem[5:cut].strip()
    #     self.seq_3_5_cuted = map(cut_3_5_adaptors, [elem for elem in self.input_file_readlines
    #                                                 if filtr_3_5_adaptors(elem, adaptor3, adaptor5)])
    #     self.seq_3_5_cuted = [elem for elem in self.seq_3_5_cuted if elem]


    def dividing_data_for_clustering(self):

        self.seq_3_5_cuted = self.input_file_readlines

        length = len(self.seq_3_5_cuted)
        a_list = range(0, length, length / self.number_of_parts)

        size_list = [[a_list[elem], a_list[elem + 1]] for elem in xrange(len(a_list) - 1)]
        size_list[-1][-1] = length

        data_for_clustering = [self.seq_3_5_cuted[start: end] for start, end in size_list]

        return data_for_clustering

    def data_clustering(self, data):
        done = []
        part = []

        def copy_seq():
            part.append([data.count(seq), seq])
            done.append(seq)

        [copy_seq() for seq in data if seq not in done]
        self.clustered_data.append(part)

    @staticmethod
    def combine_two_clustered_data(output, fragment1, fragment2):

        part1 = set(elem[1] for elem in fragment1)
        part2 = set(elem[1] for elem in fragment2)

        common = part1.intersection(part1, part2)
        not_common1 = part1 - part2
        not_common2 = part2 - part1

        part = []

        part.extend([[num, seq] for num, seq in fragment1 if seq in not_common1])
        part.extend([[num, seq] for num, seq in fragment2 if seq in not_common2])

        part.extend([[num1 + num2, seq1]
                     for num1, seq1 in ((num, seq) for num, seq in fragment1 if seq in common)
                     for num2, seq2 in ((num, seq) for num, seq in fragment2 if seq in common)
                     if seq1 == seq2])

        output.append(part)

    def combine_all_clustered_data(self):

        size_list = ((elem, elem + 1) for elem in xrange(0, len(self.clustered_data), 2))
        for x, y in size_list:
            self.combine_two_clustered_data(self.clustered_collected_A,
                                            self.clustered_data[x],
                                            self.clustered_data[y])

        size_list = ((elem, elem + 1) for elem in xrange(0, len(self.clustered_collected_A), 2))
        for x, y in size_list:
            self.combine_two_clustered_data(self.clustered_collected_B,
                                            self.clustered_collected_A[x],
                                            self.clustered_collected_A[y])

        size_list = ((elem, elem + 1) for elem in xrange(0, len(self.clustered_collected_B), 2))
        for x, y in size_list:
            self.combine_two_clustered_data(self.clustered_collected_C,
                                            self.clustered_collected_B[x],
                                            self.clustered_collected_B[y])

        print "spoko... nie zasnąłem, jestem gdzieś w połowie składania ..."
        self.combine_two_clustered_data(self.clustered_collected_D,
                                        self.clustered_collected_C[0],
                                        self.clustered_collected_C[1])

        self.clustered_collected_D = self.clustered_collected_D[0]

    def data_sorting(self, sorting):

        if sorting == "L":
            self.clustered_collected_D.sort(key=lambda x: int(x[0]), reverse=True)
        if sorting == "S":
            self.clustered_collected_D.sort(key=lambda x: len(x[1]), reverse=True)

    def save_file(self, output_file_path):

        output_string = "%s\n" \
                        "liczba indywidualnych sekwencji --> %d\n" \
                        "liczba sekwencji po zliczeniu --> %d\n\n" \
                        "Lp\tseq count\tseq length\tseq\n" % (output_file_path,
                                                            len(self.seq_3_5_cuted),
                                                            len(self.clustered_collected_D))

        for n, (num, seq) in enumerate(self.clustered_collected_D):
            output_string = "%s%d\t%d\t%d\t%s\n" % (output_string, n+1, num, len(seq), seq)

        output_file = open(output_file_path, "w")
        output_file.write(output_string)
        output_file.close()

if __name__ == '__main__':

    intro = u"Witam w programie SimpleDegradomeSequenceAnalyser.\n" \
            "Program przeprowadza wstępną analizę danych degradomowych z plików tekstowych.\n".encode(encoding="UTF-8")

    print intro
    stan = True

    while stan:
        user_input_file = raw_input("wprowadź ścieżkę do pliku z danymi input --> ")
        user_output_file = raw_input("wprowadź ścieżkę do pliku do zapisu danych output --> ")
        user_adaptor_3 = raw_input("wprowadź sekwencję adaptora 3' --> ")
        user_adaptor_5 = raw_input("wprowadź sekwencję adaptora 5' --> ")

        print "Wybierz sposób sortowania danych: \n" \
              "sortowanie danych po liczbie zliczeń sekwencji -> wpisz L\n" \
              "sortowanie danych po długości sekwencji -> wpisz S"

        sortowanie = raw_input("Sposób sortowania --> ")
        while sortowanie not in "LS":

            print "Wybierz sposób sortowania danych: \n" \
                  "sortowanie danych po liczbie zliczeń sekwencji -> wpisz L\n" \
                  "sortowanie danych po długości sekwencji -> wpisz S"
            sortowanie = raw_input("Sposób sortowania --> ")

        print "Zatem do dzieła ..."

        start = datetime.datetime.now()
        print start

        print "inicjalizacja ..."
        program = DegradomeSequenceAnalyser()

        print "otwieranie pliku input ... ", datetime.datetime.now() - start
        program.open_input_file(user_input_file)

        print "usuwanie adaptorów ... ", datetime.datetime.now() - start
        program.filtr_cut_3_5_adaptors(user_adaptor_3, user_adaptor_5)   # to trzeba zmianić

        print "dzielenie danych na 16 fragmentów oraz zliczanie sekwencji ... ", datetime.datetime.now() - start
        for num, elem in enumerate(program.dividing_data_for_clustering()):
            program.data_clustering(elem)

        print "składanie danych w całość ... to może trochę zająć ... ", datetime.datetime.now() - start
        program.combine_all_clustered_data()

        print "sortowanie danych ... ", datetime.datetime.now() - start
        program.data_sorting(sortowanie)

        print "tworzenie danych do zapisu oraz zapis danych do pliku ... ", datetime.datetime.now() - start
        program.save_file(user_output_file)

        print "WORK COMPLETED ;-)", datetime.datetime.now() - start

        print "Co dalej ? Analizujemy kolejny plik ?"

        what_next = raw_input("Wpisz T lub N: ")
        while what_next not in "TN":
            print "Co dalej ? Analizujemy kolejny plik ?"
            what_next = raw_input("Wpisz T jeśli chcesz analizować kolejny plik lub \n"
                                  "N jeśli wychodzisz z programu: ")
        if what_next == "N":
            stan = False

    exit()

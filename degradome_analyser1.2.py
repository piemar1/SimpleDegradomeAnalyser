# ! /usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import sys

__author__ = 'Marcin Pieczyński'

"""
Program służący do wstępnej analizy danych z sekwencjonowania degradomów.
Jako plik input przyjmuje plik tekstowe będące wynikiem sekwencjonowania.
Jako plik output zwraca plik z przefiltrowanymi, policzonymi i posortowanymi
sekwencjami według zadanym parametrów.

Zmiany w kodzie mogą wpłynąć na jego funkcjonowanie.
Korzystaj rozważnie.
Autor programu - Marcin Pieczyński, piemar@amu.edu.pl
"""


class DegradomeSequenceAnalyser(object):
    """
    Klasa zawierająca wszystkie metody programu potrzebne do działania programu.
    """
    def __init__(self):
        """
        Metoda inicjalizująca klasę DegradomeSequenceAnalyser.
        """
        self.input_file_readlines = []
        self.counted_data = []

        self.AA, self.AG, self.AT, self.AC = [], [], [], []
        self.GG, self.GA, self.GT, self.GC = [], [], [], []
        self.CC, self.CA, self.CG, self.CT = [], [], [], []
        self.TT, self.TG, self.TC, self.TA = [], [], [], []

        self.all_seq = [self.AA, self.AG, self.AT, self.AC,
                        self.GG, self.GA, self.GT, self.GC,
                        self.CC, self.CA, self.CG, self.CT,
                        self.TT, self.TG, self.TC, self.TA]

        self.output_file = None

    def open_input_file(self, input_file_path):
        """
        Metoda otwierająca plik input z danymi.
        """
        f = open(input_file_path)
        self.input_file_readlines = f.readlines()
        f.close()

    @staticmethod
    def filtr_adaptors(seq, adaptor):
        """
        Metoda sprawdzająca czy w zadanej sekwencji znajduje sie zadany adaptor.
        """
        if adaptor in seq:
            return True
        return False

    def filtr_cut_5_adaptors(self, adaptor):
        """
        Metoda służy do odcięcie z sekwencji adaptora3' oraz fragmentu znajdującego sie za adaptorem.
        """

        def cut_5_adaptors(seq):
            """
            Metoda bezposrednio dokonujaca odciecia adaptora.
            """
            cut = seq.index(adaptor)
            return seq[cut + len(adaptor):].strip()

        seq_5_cuted = map(cut_5_adaptors,
                          [elem for elem in self.input_file_readlines if self.filtr_adaptors(elem, adaptor)])

        self.input_file_readlines = [elem for elem in seq_5_cuted if elem]

    def filtr_cut_3_adaptors(self, adaptor):
        """
        Metoda sluzy do odciecie z sekwencji adaptora5' oraz fragmentu znajdujacego sie za adaptorem.
        """

        def cut_3_adaptors(seq):
            """
            Metoda bezpośrednio dokonująca odcięcia adaptora.
            """
            cut = seq.index(adaptor)
            return seq[:cut].strip()

        seq_3_cuted = map(cut_3_adaptors,
                          [elem for elem in self.input_file_readlines if self.filtr_adaptors(elem, adaptor)])

        self.input_file_readlines = [elem for elem in seq_3_cuted if elem]

    def dividing_data(self):
        """
        Metoda dokonująca podziału sekwencji na 16 zbiorów,
        w których każda sekwencja zaczynająca sie od innego dinukleotydu.
        """
        for seq in self.input_file_readlines:
            short = seq[:2]
            if short == "AA":
                self.AA.append(seq)
            elif short == "AG":
                self.AG.append(seq)
            elif short == "AT":
                self.AT.append(seq)
            elif short == "AC":
                self.AC.append(seq)

            elif short == "GG":
                self.GG.append(seq)
            elif short == "GA":
                self.GA.append(seq)
            elif short == "GT":
                self.GT.append(seq)
            elif short == "GC":
                self.GC.append(seq)

            elif short == "CC":
                self.CC.append(seq)
            elif short == "CA":
                self.CA.append(seq)
            elif short == "CG":
                self.CG.append(seq)
            elif short == "CT":
                self.CT.append(seq)

            elif short == "TT":
                self.TT.append(seq)
            elif short == "TG":
                self.TG.append(seq)
            elif short == "TC":
                self.TC.append(seq)
            elif short == "TA":
                self.TA.append(seq)

    def seq_counting(self, data):
        """
        Metoda dokonująca zliczeń liczby wystąpień dla poszczególnych sekwencji.
        """
        done = set()

        def copy_seq():
            self.counted_data.append([data.count(seq), seq])
            done.add(seq)

        [copy_seq() for seq in data if seq not in done]

    def data_sorting(self, sorting):
        """
        Metoda dokonująca sortowania danych sekwencji według zadanych parametrów:
        liczby zliczeń "L" lub długości sekwencji "S.
        """
        if sorting == "L":
            self.counted_data.sort(key=lambda x: int(x[0]), reverse=True)
        if sorting == "S":
            self.counted_data.sort(key=lambda x: len(x[1]), reverse=True)

    def create_output_file(self, output_file_path):
        """
        Metoda tworząca plik do zapisu danych output.
        """
        self.output_file = open(output_file_path, "w")

    def save_file(self, input_file_path, adaptor5, adaptor3):
        """
        Metoda przygotowująca dane do zapisu oraz dokonująca zapisu danych output w pliku.
        """
        if not adaptor5:
            adaptor5 = "brak"

        if not adaptor3:
            adaptor3 = "brak"

        output_string = "plik źródłowy --> {}\n" \
                        "liczba indywidualnych sekwencji --> {}\n" \
                        "liczba sekwencji po zliczeniu --> {}\n" \
                        "adaptor 5' --> {}\n" \
                        "adaptor 3' --> {}\n\n" \
                        "Lp\tseq count\tseq length\tseq\n".format(input_file_path,
                                                                  len(self.input_file_readlines),
                                                                  len(self.counted_data),
                                                                  adaptor5,
                                                                  adaptor3)

        a_string = "".join(["".join([str(no+1), "\t", str(num), "\t", str(len(seq)), "\t", seq, "\n"])
                          for no, (num, seq) in enumerate(self.counted_data)])

        output_string = output_string + a_string

        self.output_file.write(output_string)
        self.output_file.close()


if __name__ == '__main__':

    # Poniżej znajduje się zapis interfejsu tekstowego dla użytkowników programu.

    def pl(txt):
        return txt.decode(encoding='utf-8')

    def pl2(txt):
        return txt.encode(encoding=sys.stdout.encoding)

    print pl("\nWitam w programie SimpleDegradomeSequenceAnalyser 1.2.\n"
             "Program przeprowadza wstępną analizę danych degradomowych z plików tekstowych.")

    # właściwa pętla programu
    stan = True
    while stan:

        program = DegradomeSequenceAnalyser()
        print pl("\ninicjalizacja programu ...")

        # otwieranie pliku input
        file = False
        while not file:
            user_input_file = raw_input(pl2(u"\nwprowadź ścieżkę do pliku z danymi input --> ")).strip()
            print pl("\notwieranie pliku input ... ")
            try:
                program.open_input_file(user_input_file)
                file = True
            except IOError:
                print pl("\nSłuchaj no, z tym plikiem input jest coś nie tak.\n"
                         "1 - sprawdź ścieżkę dostępu do pliku\n"
                         "2 - upewnij się że plik, który chcesz otworzyć jest plikiem "
                         "tekstowym o rozwinięciu np: txt, fa, fasta, fastq")

        # otwieranie pliku output
        file = False
        while not file:
            user_output_file = raw_input(pl2(u"\nwprowadź ścieżkę do pliku dla zapisu danych output --> ")).strip()
            print pl("\ninicjalizacja pliku output ...")
            try:
                program.create_output_file(user_output_file)
                file = True
            except IOError:
                print pl("Słuchaj no, z tym plikiem output jest coś nie tak.\n"
                         "Sprawdź ścieżkę dostępu do pliku.")

        # wprowadzanie sekwencji adapterów
        user_adaptor_3 = raw_input(pl2(u"\nwprowadź sekwencję adaptora 3' --> ")).strip().upper()
        user_adaptor_5 = raw_input(pl2(u"\nwprowadź sekwencję adaptora 5' --> ")).strip().upper()

        while not user_adaptor_3 and not user_adaptor_5:
            print pl("\nProgram wymaga wprowadzenia sekwencji przynajmniej jednago z adapterów.")
            user_adaptor_3 = raw_input(pl2(u"\nwprowadź sekwencję adaptora 3' --> ")).strip().upper()
            user_adaptor_5 = raw_input(pl2(u"\nwprowadź sekwencję adaptora 5' --> ")).strip().upper()

        # wybór sposobu sortowania danych w pliku wyjściowym
        print pl("\nWybierz sposób sortowania danych: \n"
                 "sortowanie danych po liczbie zliczeń sekwencji -> wpisz L\n"
                 "sortowanie danych po długości sekwencji -> wpisz S")
        sorting_way = raw_input(pl2(u"Sposób sortowania --> ")).strip()

        while sorting_way not in "LS" or not sorting_way:
            print pl("\nWybierz sposób sortowania danych: \n"
                     "sortowanie danych po liczbie zliczeń sekwencji -> wpisz L\n"
                     "sortowanie danych po długości sekwencji -> wpisz S")
            sorting_way = raw_input(pl2(u"Sposób sortowania --> "))

        # Właściwa analiza danych !!!
        print pl("\nWszystkie dane zebrane, zatem do dzieła ... \n")
        start = datetime.datetime.now()
        print start

        if user_adaptor_3:
            print pl("usuwanie adaptorów 3' ... "), datetime.datetime.now() - start
            program.filtr_cut_3_adaptors(user_adaptor_3)

        if user_adaptor_5:
            print pl("usuwanie adaptorów 5' ... "), datetime.datetime.now() - start
            program.filtr_cut_5_adaptors(user_adaptor_5)

        print pl("dzielenie danych na 16 fragmentów ... "), datetime.datetime.now() - start
        program.dividing_data()

        print pl("liczenie sekwencji ... "), datetime.datetime.now() - start
        
        #################
        # for elem in program.all_seq:
        #     print "len(elem)", len(elem)
        #################
        
        for n, elem in enumerate(program.all_seq):
            print pl("zliczanie sekwencji w fragmencie nr - "), n + 1, "\t", datetime.datetime.now() - start
            program.seq_counting(elem)

        print pl("sortowanie sekwencji"), datetime.datetime.now() - start
        program.data_sorting(sorting_way)

        print pl("tworzenie danych do zapisu oraz zapis danych do pliku ... "), datetime.datetime.now() - start
        program.save_file(user_input_file, user_adaptor_5, user_adaptor_3)

        print pl("\nWORK COMPLETED ;-)"), datetime.datetime.now() - start

        # co dalej po zakończeniu analizy danych ?
        print pl("\nCo dalej ? Analizujemy kolejny plik ?")
        what_next = raw_input(pl2(u"Wpisz T lub N: "))
        while what_next not in "TN":

            print pl("\nHej zdecyduj się, analizujemy kolejny plik czy jednak kawa?")
            what_next = raw_input(pl2(u"Wpisz T jeśli chcesz analizować kolejny plik lub \n"
                                      u"N jeśli wychodzisz z programu: "))
        if what_next == "N":
            stan = False

    exit()

# ! /usr/bin/env python

import datetime

__author__ = 'Marcin Pieczynski'

"""
Program sluzacy do wstepnej analizy danych z sekwencjonowani degradomow.
Jako plik input przyjmuje plik tekstowe bedace wynikiem sekwencjonowania.
Jako plik output zwraca plik z przefiltrowanymi, policzonymi i posotrowanymi
sekwencjami wedlog zadanym parametrow.

Zmiany w kodzie moga wplyna na jego funkcjonowanie.
Korzystaj rozwaznie.
Autor programu - Marcin Pieczynski, piemar@amu.edu.pl
"""


class DegradomeSequenceAnalyser(object):
    """
    Klasa zawierajaca wszystkie metody programu potrzebne do dzialania.
    """
    def __init__(self):
        """
        Metoda inicjalizujaca klase DegradomeSequenceAnalyser.
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

    def open_input_file(self, input_file_path):
        """
        Metoda otwierajaca plik input z danymi.
        """
        self.input_file_readlines = open(input_file_path).readlines()

    @staticmethod
    def filtr_adaptors(seq, adaptor):
        """
        Metoda sprawdzajaca czy w zadanej sekwencji znajduje sie zadany adaptor
        """
        if adaptor in seq:
            return True
        return False

    def filtr_cut_3_adaptors(self, adaptor3):
        """
        Metoda sluzy do odciecie z sekwencji adaptora3' oraz fragmentu znajdujacego sie za adaptorem.
        """

        def cut_3_adaptors(seq):
            """
            Metoda bezposrednio dokonujaca odciecia adaptora
            """
            cut = seq.index(adaptor3)
            return seq[cut + len(adaptor3):].strip()

        seq_3_cuted = map(cut_3_adaptors,
                          [elem for elem in self.input_file_readlines if self.filtr_adaptors(elem, adaptor3)])

        self.input_file_readlines = [elem for elem in seq_3_cuted if elem]

    def filtr_cut_5_adaptors(self, adaptor5):
        """
        Metoda sluzy do odciecie z sekwencji adaptora5' oraz fragmentu znajdujacego sie za adaptorem.
        """

        def cut_5_adaptors(seq):
            """
            Metoda bezposrednio dokonujaca odciecia adaptora
            """
            cut = seq.index(adaptor5)
            return seq[:cut].strip()

        seq_5_cuted = map(cut_5_adaptors,
                          [elem for elem in self.input_file_readlines if self.filtr_adaptors(elem, adaptor5)])

        self.input_file_readlines = [elem for elem in seq_5_cuted if elem]

    def dividing_data(self):
        """
        Metoda dokonujaca podzialu sekwencji na 16 zbiorow,
        w ktorych kazda sekwencja zaczynajaca sie od innego dinukleotydu.
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
        Metoda dokonujaca zliczen liczby wystapien dla poszczegolnych sekwencji.
        """
        done = set()

        def copy_seq():
            self.counted_data.append([data.count(seq), seq])
            done.add(seq)

        [copy_seq() for seq in data if seq not in done]
        # self.counted_data.extend(part)

    def data_sorting(self, sorting):
        """
        Metoda dokonujaca sortowania danych sekwencji wedlog zadanych parametrow:
        liczby zliczen "L" lub dlugosci sekwencji "S.
        """

        if sorting == "L":
            self.counted_data.sort(key=lambda x: int(x[0]), reverse=True)
        if sorting == "S":
            self.counted_data.sort(key=lambda x: len(x[1]), reverse=True)

    def save_file(self, output_file_path):
        """
        Metoda przygotowujaca dane do zapisu oraz dokonujaca zapisu danych output w pliku.
        """

        output_string = "%s\n" \
                        "liczba indywidualnych sekwencji --> %d\n" \
                        "liczba sekwencji po zliczeniu --> %d\n\n" \
                        "Lp\tseq count\tseq length\tseq\n" % (output_file_path,
                                                              len(self.input_file_readlines),
                                                              len(self.counted_data))

        for no, (num, seq) in enumerate(self.counted_data):
            output_string = "%s%d\t%d\t%d\t%s\n" % (output_string, no+1, num, len(seq), seq)

        output_file = open(output_file_path, "w")
        output_file.write(output_string)
        output_file.close()


if __name__ == '__main__':

    # Ponizej znajduje sie zapis interfajsu tekstowego dla uzytkownikow programu.

    intro = "Witam w programie SimpleDegradomeSequenceAnalyser.\n" \
            "Program przeprowadza wstepna analize danych degradomowych z plikow tekstowych.\n"

    print intro
    stan = True

    while stan:
        user_input_file = raw_input("wprowadz sciezke do pliku z danymi input --> ").strip()
        user_output_file = raw_input("wprowadz sciezke do pliku do zapisu danych output --> ").strip()

        user_adaptor_3 = raw_input("wprowadz sekwencje adaptora 3' --> ").strip()
        user_adaptor_5 = raw_input("wprowadz sekwencje adaptora 5' --> ").strip()

        while not user_adaptor_3 and not user_adaptor_5:
            print "\nProgram wymaga wprowadzenia sekwencji przynajmniej jednago z adapterow."
            user_adaptor_3 = raw_input("wprowadz sekwencje adaptora 3' --> ").strip()
            user_adaptor_5 = raw_input("wprowadz sekwencje adaptora 5' --> ").strip()

        print "\nWybierz sposob sortowania danych: \n" \
              "sortowanie danych po liczbie zliczen sekwencji -> wpisz L\n" \
              "sortowanie danych po dlugosci sekwencji -> wpisz S"

        sorting_way = raw_input("Sposob sortowania --> ").strip()

        while sorting_way not in "LS" or not sorting_way:

            print "\nWybierz sposob sortowania danych: \n" \
                  "sortowanie danych po liczbie zliczen sekwencji -> wpisz L\n" \
                  "sortowanie danych po dlugosci sekwencji -> wpisz S"
            sorting_way = raw_input("Sposob sortowania --> ")

        print "\nZatem do dziela ... \n"
        start = datetime.datetime.now()
        print start

        print "inicjalizacja ..."
        program = DegradomeSequenceAnalyser()

        print "otwieranie pliku input ... ", datetime.datetime.now() - start
        program.open_input_file(user_input_file)

        if user_adaptor_3:
            print "usuwanie adaptorow 3' ... ", datetime.datetime.now() - start
            program.filtr_cut_3_adaptors(user_adaptor_3)

        if user_adaptor_5:
            print "usuwanie adaptorow 5' ... ", datetime.datetime.now() - start
            program.filtr_cut_5_adaptors(user_adaptor_5)

        print "dzielenie danych na 16 fragmentow ... ", datetime.datetime.now() - start
        program.dividing_data()

        print "liczenie sekwencji ... ", datetime.datetime.now() - start
        for n, elem in enumerate(program.all_seq):
            print "zliczanie sekwencji w fragmencie nr - ", n + 1
            program.seq_counting(elem)

        print "sortowanie sekwencji", datetime.datetime.now() - start
        program.data_sorting(sorting_way)

        print "tworzenie danych do zapisu oraz zapis danych do pliku ... ", datetime.datetime.now() - start
        program.save_file(user_output_file)

        print "WORK COMPLETED ;-)", datetime.datetime.now() - start

        print "Co dalej ? Analizujemy kolejny plik ?"

        what_next = raw_input("Wpisz T lub N: ")
        while what_next not in "TN":

            print "\nHej zdecyduj sie, analizujemy kolejny plik czy jednak kawa?"
            what_next = raw_input("Wpisz T jesli chcesz analizowac kolejny plik lub \n"
                                  "N jesli wychodzisz z programu: ")
        if what_next == "N":
            stan = False

    exit()

# ! /usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import sys
from string import ascii_uppercase as up

__author__ = 'Marcin Pieczyński'

"""
Program służący do wstępnej analizy danych z sekwencjonowania degradomów.
Jako plik input przyjmuje pliki uzyskane dzięki programowi degradome_analyser.
Zalecana wersja programu degradome_analyser1.2.py

Jako plik output zwraca plik z zawierqjący informacje o liczbie wystąpień sekwencji w poszczególnych plikach input.

Zmiany w kodzie mogą wpłynąć na jego funkcjonowanie.
Korzystaj rozważnie.
Autor programu - Marcin Pieczyński, piemar@amu.edu.pl
"""

class DegradomeResultsAssambler(object):
    def __init__(self):

        self.input_all_data = []
        self.spec_data = []
        self.final_data = []
        self.output_file = None
        self.size = 0

    def open_input_file(self, input_file_path):
        """
        Metoda otwierająca plik input z danymi.
        """
        f = open(input_file_path)
        self.input_all_data.append(f.readlines())
        f.close()

    def read_spec_data(self):
        """
        Metoda odczytujący dane z poszczególnych plików input.
        """
        for set in self.input_all_data:
            for n, line in enumerate(set):
                if line[:2] == "Lp":
                    number = n + 1
            set = set[number:]

            data = []
            for line in set:
                line = line.split("\t")
                data.append([line[1], line[3].strip()])

            self.spec_data.append(data)
        self.size = len(self.spec_data)

    def create_final_data(self):
        """
        Metoda inicjalizuje finalne wyniki w postaci stringa.
        """
        for elem in self.spec_data[0]:
            d = [str(elem[0]), str(elem[1])]
            for x in xrange(self.size - 1):
                d.insert(1, '0')
            self.final_data.append(d)

    def data_assambling(self, spec_data_number):
        """
        Metoda uzupełnia finalne wyniki danymi z poszczególnych plików input.
        """

        setA = set(seq for num,seq in self.spec_data[0])
        setB = set(seq for num,seq in self.spec_data[spec_data_number])

        common = setA.intersection(setA, setB)   # seq wspólne
        not_common2 = setB - setA                # seq z B jakie trezba dodać do A

        print "liczba sekwencji wspólnych dla plików    ", len(common)
        print "liczba wszystkich sekwencji w plikach    ", len(not_common2) + len(common)

        for no, (num, seq) in enumerate(self.spec_data[spec_data_number]):
            # print no
            if seq in not_common2:
                self.add_not_common_data(num, seq, spec_data_number)
                # funkcja dodająca dane na końcu listy

            elif seq in common:
                # print no
                self.add_common_data(num, seq, spec_data_number)
                # funkcja dodająca dane w środku listy

    def add_not_common_data(self, num, seq, spec_data_number):
        """
        Metoda uzupełniająca dane dla sekwencji nie wspólnych dla analizowanych plików.
        """
        data = ['0' for x in xrange(self.size)]
        data.extend([seq])

        data[spec_data_number] = num
        self.final_data.append(data)

    def add_common_data(self, num, seq, spec_data_number):
        """
        Metoda uzupełniająca dane dla sekwencji wspólnych dla analizowanych plików.
        """
        for n, line in enumerate(self.final_data):
            if seq in line:
                self.final_data[n][spec_data_number] = num
                break

    def create_output_file(self, output_file_path):
        """
        Metoda tworząca plik do zapisu danych output.
        """
        self.output_file = open(output_file_path, "w")

    def save_file(self, output_file_path, input_files_paths):
        """
        Metoda przygotowująca dane do zapisu oraz dokonująca zapisu danych output w pliku.
        """
        output_string = ""
        for number, file_path in enumerate(input_files_paths):
            output_string += "{} - \t{}\n".format(str(up[number]), file_path)

        output_string += "\n\nLp"
        for number, file_path in enumerate(input_files_paths):
            output_string = "{}\t{}".format(output_string, up[number])

        output_string += "\tseq\n"

        a_string = "\n".join(["".join([str(no+1), "\t", "\t".join(line)])
                              for no, (line) in enumerate(self.final_data)])

        output_string += a_string
        outputfile = open(output_file_path, "w")
        outputfile.write(output_string)
        outputfile.close()


if __name__ == '__main__':

    def pl(txt):
        return txt.decode(encoding='utf-8')

    def pl2(txt):
        return txt.encode(encoding=sys.stdout.encoding)

    print pl("\nWitam w programie DegradomeResultsAssambler 1.0.\n"
             "Program zbiera dane zawarte w plikach wygenerowanych za pomocą programu\n"
             "SimpleDegradomeSequenceAnalyser (preferowane wersja - 1.2).\n"
             "Jako output program zwraca plik zawierający informacje o występowaniu\n"
             "poszczególnych sekwencji w plikach input.\n")

    file_number = raw_input("Jaką liczbę plików chcesz analizować, wpisz liczbę --> ")
    number = False
    while not number:
        try:
            file_number = int(file_number)
            if file_number > 1:
                number = True
            else:
                print """Podano nieprawidłową liczbę plików do analizy.
                Prawidłowa liczba plików to 2 i więcej."""
                file_number = raw_input("Jaką liczbę plików chcesz analizować, wpisz liczbę --> ")

        except ValueError:
            print "Należy wpisać cyfrę a nie tekst."
            file_number = raw_input("Jaką liczbę plików chcesz analizować, wpisz liczbę --> ")

    program = DegradomeResultsAssambler()
    program.input_files_paths = []
    for number in xrange(file_number):
        file = False
        while not file:
            f = raw_input("Podaj ścieżke do pliku input--> ")
            print pl("otwieranie pliku input ... \n")
            try:
                program.open_input_file(f)
                program.input_files_paths.append(f)
                file = True
            except IOError:
                print pl("Słuchaj no, z tym plikiem input jest coś nie tak.\n"
                         "1 - sprawdź ścieżkę dostępu do pliku\n"
                         "2 - upewnij się że plik, który chcesz otworzyć jest plikiem "
                         "tekstowym o rozwinięciu np: txt, fa, fasta, fastq")

    # otwieranie pliku output
    file = False
    while not file:
        user_output_file = raw_input(pl2(u"Wprowadź ścieżkę do pliku dla zapisu danych output --> ")).strip()
        print pl("\ninicjalizacja pliku output ...")
        try:
            program.create_output_file(user_output_file)
            file = True
        except IOError:
            print pl("Słuchaj no, z tym plikiem output jest coś nie tak.\n"
                     "Sprawdź ścieżkę dostępu do pliku.")

    start = datetime.datetime.now()
    print "Start, godzina --> ", start

    program.read_spec_data()
    program.create_final_data()

    print datetime.datetime.now() - start

    for elem in xrange(1, program.size):
        print "Assambling", datetime.datetime.now() - start
        print "Zbieranie danych dla analizowanej pary plików..."
        program.data_assambling(elem)

    print "zapis danych do pliku ...    ", datetime.datetime.now() - start

    program.save_file(user_output_file,  program.input_files_paths)

    print"Work Completed ;-)"
    print datetime.datetime.now() - start

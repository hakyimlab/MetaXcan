__author__ = 'heroico'

import csv

class STFS(object):
    """sample file table format"""
    ID = 0
    POP = 1
    GROUP = 2
    SEX = 3

class Person(object):
    """A person."""
    def __init__(self):
        self.id = None
        self.population = None
        self.group = None
        self.sex = None

    def toTextLine(self):
        return " ".join([self.id, self.population, self.group, self.sex])

    @classmethod
    def loadPersonFromSampleRow(cls, row):
        person = Person()
        person.id = row[STFS.ID]
        person.population = row[STFS.POP]
        person.group = row[STFS.GROUP]
        person.sex = row[STFS.SEX]
        return  person

    @classmethod
    def loadPersonFromSampleRowIfFilter(cls, row, population_filters, individual_filters):
        person = None

        if not "None" in population_filters and\
            row[STFS.GROUP] not in population_filters:
            return person

        if len(individual_filters):
            id = row[STFS.ID]
            match = None
            for filter in individual_filters:
                match = filter.match(id)
                if match:
                    break
            if not match:
                return person

        person = Person.loadPersonFromSampleRow(row)
        return person

    @classmethod
    def loadPeople(cls, input, delim=' ', skip_header=True):
        people = []
        with open(input, 'rb') as csv_file:
            reader = csv.reader(csv_file, delimiter=delim, quotechar='"')
            for row in reader:
                if skip_header and reader.line_num == 1:
                    continue
                person = Person.loadPersonFromSampleRow(row)
                people.append(person)
        return people

    @classmethod
    def loadFilteredPeople(cls, input_path, population_filters = ["EUR"], individual_filters =[], row_delimiter=' ', skip_header=True):
        filtered = []
        with open(input_path, 'rb') as csv_file:
            reader = csv.reader(csv_file, delimiter=row_delimiter, quotechar='"')
            for row in reader:
                if skip_header and reader.line_num == 1:
                    continue

                person = Person.loadPersonFromSampleRowIfFilter(row, population_filters, individual_filters)
                if person is not None:
                    filtered.append(person)
        return filtered

    @classmethod
    def buildFilteredSamples(cls, input_path, output, population_filters = ["EUR"], individual_filters =[], row_delimiter=' ', skip_header=True):
        filtered = cls.loadFilteredPeople(input_path, population_filters, individual_filters, row_delimiter, skip_header)

        with open(output, 'w+') as output_file:
            output_file.write(" ".join(["ID", "POP", "GROUP", "SEX"])+"\n")
            for person in filtered:
                output_file.write(person.toTextLine()+"\n")

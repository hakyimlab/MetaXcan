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
    def __init__(self, id=None, population=None, group=None, sex=None):
        self.id = id
        self.population = population
        self.group = group
        self.sex = sex

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
    def loadPersonFromSampleRowIfFilter(cls, row, group_filters, individual_filters):
        person = None

        if not "None" in group_filters and\
            row[STFS.GROUP] not in group_filters:
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
        with open(input, "r") as csv_file:
            reader = csv.reader(csv_file, delimiter=delim, quotechar='"')
            for row in reader:
                if skip_header and reader.line_num == 1:
                    continue
                person = Person.loadPersonFromSampleRow(row)
                people.append(person)
        return people

    @classmethod
    def loadFilteredPeople(cls, input_path, group_filters = ["EUR"], individual_filters =[], row_delimiter=' ', skip_header=True):
        filtered = []
        with open(input_path, 'r') as csv_file:
            reader = csv.reader(csv_file, delimiter=row_delimiter, quotechar='"')
            for row in reader:
                if skip_header and reader.line_num == 1:
                    continue

                person = Person.loadPersonFromSampleRowIfFilter(row, group_filters, individual_filters)
                if person is not None:
                    filtered.append(person)
        return filtered

    @classmethod
    def buildFilteredSamples(cls, input_path, output, group_filters = ["EUR"], individual_filters =[], row_delimiter=' ', skip_header=True):
        filtered = cls.loadFilteredPeople(input_path, group_filters, individual_filters, row_delimiter, skip_header)

        with open(output, 'w+') as output_file:
            output_file.write(" ".join(["ID", "POP", "GROUP", "SEX"])+"\n")
            for person in filtered:
                output_file.write(person.toTextLine()+"\n")

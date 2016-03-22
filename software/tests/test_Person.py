import unittest
import sys
import os
import shutil
import re

if "DEBUG" in sys.argv:
    sys.path.insert(0, "..")
    sys.path.insert(0, "../../")
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

from metax.Person import Person

class TestPerson(unittest.TestCase):
    def testPersonConstructor(self):
        p = Person("a", "b", "c", "d")
        self.assertPerson(p)

    def testPersonToLine(self):
        p = Person("a", "b", "c", "d")
        l = p.toTextLine()
        self.assertEqual("a b c d", l)

    def testPersonLoadFromRow(self):
        row = ["a", "b", "c", "d"]
        p = Person.loadPersonFromSampleRow(row)
        self.assertPerson(p)

    def testPersonLoadFromRowIfFilter(self):
        row1 = ["a", "b", "c", "d"]
        p = Person.loadPersonFromSampleRowIfFilter(row1, group_filters=["None", "EUR"], individual_filters=[])
        self.assertPerson(p)

        p = Person.loadPersonFromSampleRowIfFilter(row1, group_filters=["EUR"], individual_filters=[])
        self.assertIsNone(p)

        p = Person.loadPersonFromSampleRowIfFilter(row1, group_filters=["c", "EUR"], individual_filters=[])
        self.assertPerson(p)

        p = Person.loadPersonFromSampleRowIfFilter(row1, group_filters=["c", "EUR"], individual_filters=[re.compile("Nope")])
        self.assertIsNone(p)

        p = Person.loadPersonFromSampleRowIfFilter(row1, group_filters=["c", "EUR"], individual_filters=[re.compile("Nope"), re.compile("a")])
        self.assertPerson(p)

    def testLoadPeople(self):
        #Shouldn't work when wrong delimiter is specified
        thrown = False
        try:
            person = Person.loadPeople("tests/_td/dosage_set_1/set.sample", delim=",")
        except:
            thrown = True
        self.assertEqual(thrown, True)

        #default load
        people = Person.loadPeople("tests/_td/dosage_set_1/set.sample")
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 5)
        self.assertPerson(people[0], "ID1", "K", "HERO", "male")
        self.assertPerson(people[1], "ID2", "K", "HERO", "female")
        self.assertPerson(people[2], "DI5", "K", "HERO", "male")
        self.assertPerson(people[3], "ID3", "K", "HERO", "female")
        self.assertPerson(people[4], "B1", "L", "T", "female")

        #let's pretend the file header is a person entry
        people = Person.loadPeople("tests/_td/dosage_set_1/set.sample", skip_header=False)
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 6)
        self.assertPerson(people[0], "ID", "POP", "GROUP", "SEX")
        self.assertPerson(people[1], "ID1", "K", "HERO", "male")
        self.assertPerson(people[2], "ID2", "K", "HERO", "female")
        self.assertPerson(people[3], "DI5", "K", "HERO", "male")
        self.assertPerson(people[4], "ID3", "K", "HERO", "female")
        self.assertPerson(people[5], "B1", "L", "T", "female")

    def testLoadFilteredPeople(self):
        #Shouldn't work when wrong delimiter is specified
        thrown = False
        try:
            person = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", delim=",")
        except:
            thrown = True
        self.assertEqual(thrown, True)

        #default load
        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample")
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 0)

        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", group_filters=["HERO"])
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 4)
        self.assertPerson(people[0], "ID1", "K", "HERO", "male")
        self.assertPerson(people[1], "ID2", "K", "HERO", "female")
        self.assertPerson(people[2], "DI5", "K", "HERO", "male")
        self.assertPerson(people[3], "ID3", "K", "HERO", "female")

        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", group_filters=["HERO", "T"])
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 5)
        self.assertPerson(people[0], "ID1", "K", "HERO", "male")
        self.assertPerson(people[1], "ID2", "K", "HERO", "female")
        self.assertPerson(people[2], "DI5", "K", "HERO", "male")
        self.assertPerson(people[3], "ID3", "K", "HERO", "female")
        self.assertPerson(people[4], "B1", "L", "T", "female")

        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", group_filters=["T"])
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 1)
        self.assertPerson(people[0], "B1", "L", "T", "female")

        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", group_filters=["HERO"], individual_filters=[re.compile("ID")])
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 3)
        self.assertPerson(people[0], "ID1", "K", "HERO", "male")
        self.assertPerson(people[1], "ID2", "K", "HERO", "female")
        self.assertPerson(people[2], "ID3", "K", "HERO", "female")

        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", group_filters=["HERO"], individual_filters=[re.compile("ID"), re.compile("DI")])
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 4)
        self.assertPerson(people[0], "ID1", "K", "HERO", "male")
        self.assertPerson(people[1], "ID2", "K", "HERO", "female")
        self.assertPerson(people[2], "DI5", "K", "HERO", "male")
        self.assertPerson(people[3], "ID3", "K", "HERO", "female")

        people = Person.loadFilteredPeople("tests/_td/dosage_set_1/set.sample", group_filters=["HERO"], individual_filters=[re.compile("DI")])
        self.assertIsNotNone(people)
        self.assertEqual(len(people), 1)
        self.assertPerson(people[0], "DI5", "K", "HERO", "male")

    def testBuildFilteredSamples(self):
        if os.path.exists("_test"):
            shutil.rmtree("_test")
        os.makedirs("_test")

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample")
        expected = ["ID POP GROUP SEX"]
        self.assertSamplesFile("_test/s.sample", expected)

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample", group_filters=["HERO"])
        expected = ["ID POP GROUP SEX",
                    "ID1 K HERO male",
                    "ID2 K HERO female",
                    "DI5 K HERO male",
                    "ID3 K HERO female"]
        self.assertSamplesFile("_test/s.sample", expected)

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample", group_filters=["HERO", "T"])
        expected = ["ID POP GROUP SEX",
                    "ID1 K HERO male",
                    "ID2 K HERO female",
                    "DI5 K HERO male",
                    "ID3 K HERO female",
                    "B1 L T female"]
        self.assertSamplesFile("_test/s.sample", expected)

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample", group_filters=["T"])
        expected = ["ID POP GROUP SEX",
                    "B1 L T female"]
        self.assertSamplesFile("_test/s.sample", expected)

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample", group_filters=["HERO"], individual_filters=[re.compile("ID")])
        expected = ["ID POP GROUP SEX",
                    "ID1 K HERO male",
                    "ID2 K HERO female",
                    "ID3 K HERO female"]
        self.assertSamplesFile("_test/s.sample", expected)

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample", group_filters=["HERO"], individual_filters=[re.compile("DI")])
        expected = ["ID POP GROUP SEX",
                    "DI5 K HERO male"]
        self.assertSamplesFile("_test/s.sample", expected)

        Person.buildFilteredSamples("tests/_td/dosage_set_1/set.sample", "_test/s.sample", group_filters=["HERO"], individual_filters=[re.compile("DI"), re.compile("ID")])
        expected = ["ID POP GROUP SEX",
                    "ID1 K HERO male",
                    "ID2 K HERO female",
                    "DI5 K HERO male",
                    "ID3 K HERO female"]
        self.assertSamplesFile("_test/s.sample", expected)

        shutil.rmtree("_test")


    def assertPerson(self, p, id="a", population="b", group="c", sex="d"):
        self.assertIsNotNone(p)
        self.assertEqual(id, p.id)
        self.assertEqual(population, p.population)
        self.assertEqual(group, p.group)
        self.assertEqual(sex, p.sex)

    def assertSamplesFile(self, path, expected):
        with open(path) as file:
            for i, line in enumerate(file):
                e = expected[i]
                self.assertEqual(line.strip(), e)

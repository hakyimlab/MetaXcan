#!/usr/bin/env python
import pandas
import os
import re
from sqlalchemy import create_engine


engine = create_engine('sqlite:///tests/_td/dbs_2/model_sim_1.db')
weights_1 = pandas.read_sql_table('weights', engine)
weights_1.to_csv("tests/_td/dbs_2/model_sim_1.weights", index=False, sep="\t")

engine = create_engine('sqlite:///tests/_td/dbs_2/model_sim_2.db')
weights_2 = pandas.read_sql_table('weights', engine)
weights_2.to_csv("tests/_td/dbs_2/model_sim_2.weights", index=False, sep="\t")


def get_model_weights(path):
    engine = create_engine('sqlite:///'+path)
    return pandas.read_sql_table('weights', engine)

def get_weights_in_models(path):
    f = [x for x in os.listdir(path) if ".weights" in x]
    r = []
    for f_ in f:
        p_ = os.path.join(path,f_)
        weights = pandas.read_table(p_)
        weights["model"] = f_.split(".")[0]
        r.append(weights)
    return pandas.concat(r)

def models_to_weights(ip, op):
    f = os.listdir(ip)
    for f_ in f:
        p = os.path.join(ip, f_)
        w = get_model_weights(p)
        p_ = os.path.join(op, f_.split(".db")[0]+ ".weights")
        w.to_csv(p_, index=False, sep="\t")

models_to_weights("tests/_td/dbs_3", "tests/_td/dbs_3")
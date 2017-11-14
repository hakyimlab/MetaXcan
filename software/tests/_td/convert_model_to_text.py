import pandas
from sqlalchemy import create_engine

engine = create_engine('sqlite:///tests/_td/dbs_2/model_sim_1.db')
weights_1 = pandas.read_sql_table('weights', engine)
weights_1.to_csv("tests/_td/dbs_2/model_sim_1.weights", index=False, sep="\t")

engine = create_engine('sqlite:///tests/_td/dbs_2/model_sim_2.db')
weights_2 = pandas.read_sql_table('weights', engine)
weights_2.to_csv("tests/_td/dbs_2/model_sim_2.weights", index=False, sep="\t")

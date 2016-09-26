# -------------------------- EXPORTING DATA -----------------------------
from nested_dict import nested_dict
import csv
def myprint(d):
    with open('simulation_summary.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for k, v in d.iteritems():
            if isinstance(v, dict):
                myprint(v)
            else:
                writer.writerows([ k, v ])
                # print('{0} = {1}, '.format(k, v))


def myprint2(d):
    with open('simulation_summary.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for k, v in d.iteritems_flat():
            if isinstance(v, dict):
                myprint(v)
            else:
                writer.writerows([ k, v ])
                # print('{0} = {1}, '.format(k, v))

def myprint3(d):
    with open('simulation_summary.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for k, v in d.items_flat():
            if isinstance(v, dict):
                myprint(v)
            else:
                writer.writerows([ k, v ])
                # print('{0} = {1}, '.format(k, v))


from py_mentat import *

def main():
    initialCallusYoungsMod = 3
    initialCallusPoissons = 0.3

    print("NEW RUN")

    m = py_get_int("nsets()")
    callusSetID = 0
    callusElementIDs = []

    for i in range(1,m+1):
        id = py_get_int("set_id(%d)" % i)
        sn = py_get_string("set_name(%d)" % id)
        st = py_get_string("set_type(%d)" % id)
        n = py_get_int("nset_entries(%d)" % id)

        print (id,"Set ",sn,"is a ",st," set with ",n,"entries")

    for i in range(1,m+1):
        id = py_get_int("set_id(%d)" % i)
        if(py_get_string("set_name(%d)" % id) == "callus" and py_get_string("set_type(%d)" % id) == "element"):
            callusSetID = id

    numCallusElements = py_get_int("nset_entries(%d)" % callusSetID)
    for i in range(1,numCallusElements+1):
        k = py_get_int("set_entry(%d,%d)" % (callusSetID,i))
        callusElementIDs.append(k)
        print(k)

    count = 1
    for i in callusElementIDs:
        py_send("*new_mater standard *mater_option general:state:solid *mater_option general:skip_structural:off *mater_name callus%d *mater_param structural:youngs_modulus %f *mater_param structural:poissons_ratio %f *add_mater_elements %d #" % (i,initialCallusYoungsMod,initialCallusPoissons,i))

        print(count)
        count += 1

    py_send("*save_model")

if __name__ == '__main__':
    py_connect("", 40007)
    main()
    py_disconnect()


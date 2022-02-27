from asyncore import read
import json
import re
import os
import glob
import pip


class Jeison:
    def __init__(self, file):
        self.file = file
        jsonfile = open(file)
        self.data = json.load(jsonfile)
        self.service_list = list(self.data.keys())
        jsonfile.close()

    def servlist(self):
        """
        Lists all available services found in the .json file
        """
        print("Available services:")
        return self.service_list

    def servinfo(self, service):
        if service in self.service_list:
            return self.data[service]
        else:
            return "service not found"
       
    def servstr(self, service):
        """
        Returns all the data found for the given service as strings
        """
        print(service, "info:")
        servdat = self.data[service]
        for key,value in servdat.items():
            if isinstance(value, dict):
                for key2,value2 in value.items():
                    print(key+"/"+key2, ":", ', '.join(map(str, value2)))
            elif isinstance(value, list):
                print(key, ":", ', '.join(map(str, value)))
            else:
                print(key, ":", value)

    def finder(self, service, found):
        if found in self.data[service]:
            return self.data[service][found]
        else:
            for key,value in self.data[service].items():
                if isinstance(value, dict):
                    if found in self.data[service][key]:
                        return self.data[service][key][found]
            return None
        
    def get_find_aux(self, service, found):
        if isinstance(service, dict):
            if found in service:
                return service[found]
            for key in service:
                finder = self.get_find_aux(service[key], found)
                if finder is not None:
                    return finder
        if isinstance(service, list):
            for object in service:
                finder = self.get_find_aux(object, found)
                if finder is not None:
                    return finder
        return None

    def get_find_deep(self, service, found):
        if service in self.service_list:
            serv_dict = self.data[service]
            call_aux = self.get_find_aux(serv_dict, found)
            return call_aux
        else:
            return "service not found"

    def help(self):
        return f"Servlist: List the available services"

"""

            UNUSED CODE, ONLY NEEDED FOR TESTING

servicios = Jeison("/home/bioinfoadm/Desktop/Python_testing/services.json")
servicios.servstr("mtbseq")
servicios.servinfo("mtbseq")
servicios.servlist()
servicios.finder("mtbseq", "folders")


            POSIBLE FUNCION RECURSIVA

def dictfinder(nestedict, service):
    data = []
    for key,value in nestedict.items():
        if key == service:
            data.append(value)
        elif isinstance(value, dict):
            values = dictfinder(value, service)
            for value in values:
                data.append(value)
        elif isinstance(value, list):
            for object in value:
                if isinstance(object, dict):
                    objects = dictfinder(object, service)
                    for value in objects:
                        data.append(value)

jsondir = glob.glob("services.json")

jsonfile = open("/home/bioinfoadm/Desktop/Python_testing/services.json")

data = json.load(jsonfile)

service = "assembly_annotation"
what = "description"

print(data[service][what])

for key,value in data.items():
    print(key, ":", value["description"], "+", value["clean"]["folders"][0])
    x = isinstance (value["clean"], dict)
    y = isinstance (value["label"], dict)
    print(x)
    print(y)

for key,value in data["mtbseq"].items():
    if isinstance(value, dict):
        print(value)
        #print(data["mtbseq"][key])
    else:
        print("nope")


for service,key in data.items():
    for metadt,key in data (as.str(service)):
        print(metadt, ":", key)
        print("---")

    


product = jsonObject['assembly_annotation']

print(product["template"])

randomjson =  '{ "name":"John", "age":30, "city":"New York"}'

y = json.loads(randomjson)

"""

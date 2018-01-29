from JsonFile import JsonFile
import os

test = JsonFile(
    '/home/jakobh/Dokumente/PEDIA/PEDIA-workflow-new_download/1_qualityCheck/json/mapped/16057.json')
d=test.raw
print(d)
test.dumptofile('/home/jakobh/Dokumente/PEDIA/PEDIA-workflow-new_download/1_qualityCheck/json/mapped/16057.json')


def common_elements(list1, list2):
    return [element for element in list1 if element in list2]


ngap_best = 0 
glimit_best = 0
fit_best = -1

for ngap in range(5,50):
    for glimit in range(100,1500,100):
        data = t2
        data = updateArray(data)
        #data = findClusterRedshift(data)
        data['CLUSZ'] = tZ
        data = findSeperationSpatial(data, center)
        data = findLOSV(data)
        # make initial cuts
        mask = abs(data['LOSV']) < 5000
        data = data[mask]
        while True:
            try:
                if size == data.size:
                    break
            except NameError:
                pass

            size = data.size
            #print 'size', data.size

            #data = rejectInterlopers(data)
            flag = False
            try:
                x = shifty_gapper(data['SEP'], data['Z'], tZ, ngap=ngap,
                        glimit=glimit)
            except:
                flag = True
                break
            data = data[x]
            #data = findLOSVD(data)
            data = findLOSVDgmm(data)
            data['LOSVD'] = data['LOSVDgmm']

            data = findR200(data)
            mask = data['SEP'] < data['R200'][0]
            data = data[mask]

            data = findClusterRedshift(data)
            data = findSeperationSpatial(data, center)
            data = findLOSV(data)

        if not flag:
            matched = len(common_elements(t['HALOID'], data['HALOID']))
            fit = matched/t.size + 1/data.size
            if fit > fit_best:
                fit_best = fit
                ngap_best = ngap
                glimit_best = glimit
        else:
            pass



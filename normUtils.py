import statistics


def normalise(x,method = "z"):
    if method == "z":
        sd = statistics.pstdev(x)
        m = sum(x)/len(x)
        return((x - m )/sd )
    else:
        print("method " + method + " not developed yet")
    


def disc(x):
    sd = statistics.pstdev(x)
    output_series = []
    for i in x:
        if i > sd:
            output_series.append(1)
        elif i < -sd:
            output_series.append(-1)
        else:
            output_series.append(0)
    return(output_series)
    
    
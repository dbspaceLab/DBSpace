""" Return the phases and the labels associated with them"""
all_phases = (
    ["A0" + str(num) for num in range(4, 0, -1)]
    + ["B0" + str(num) for num in range(1, 5)]
    + ["C0" + str(num) for num in range(1, 10)]
    + ["C" + str(num) for num in range(10, 25)]
)


def Phase_List(exprs="all"):
    if exprs == "all":
        return all_phases
    elif exprs == "ephys":
        return all_phases[4:]
    elif exprs == "therapy":
        return all_phases[8:]
    elif exprs == "notherapy":
        return all_phases[0:8]


""" Function to plot bands since this has been annoying every time I've had to recode the thing from scratch"""

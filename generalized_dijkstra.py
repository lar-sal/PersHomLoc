# Forced path choice for sparse matrix, pick the node of lowest degree, tie-breaker is age for now
def dijkstra_solve(cols,w,t):
    s = general_dijkstra(cols, w, t, w_t=0)
    solution_size = s[1]
    solution = []
    memory_peak = "NaN"
    return solution_size, solution, "NaN", memory_peak


def face_choice(current, cols):
    node = current[0]
    minimal_neighbours = float('inf')

    for i in node:
        # Finner kolonnne som inneholder i, men som ikke har vært tidligere i stien
        neigh_cols = [(j, col) for (j, col) in enumerate(cols) if i in col and col not in current[2]]
        l = len(neigh_cols)

        # Hvis vi finner et fjes uten naboer, eller en nabo, velger vi det og slutter å lete videre
        if l == 0 or l == 1:
            next_row = i
            next_cols = neigh_cols
            break

        # Hvis det ikke er noen fjes med 0 eller 1 nabo så velger vi fjeset med færrest naboer
        elif l < minimal_neighbours:
            next_row = i
            next_cols = neigh_cols
            minimal_neighbours = l

    return (next_row, next_cols)




def cycle_in_path(path):
    """
    Check if the path already has a cycle, and can therefore not be part of a solution
    Just reduction algorithm described in Edelsbrunner-Harer which terminates if a column is zero
    Returns the zero-column index, so can be used to finding a target (column that can be summed to zero)
    """
    R= [sorted(list(p)) for p in path]
    low_one = [r[0] for r in R]
    D=[]
    while R != D:
        D = R.copy()
        for j, r in enumerate(R):
            if r == []:
                return (j)
            while low_one.index(low_one[j]) != j:
                l = low_one.index(low_one[j])
                symdif = set(R[j]).symmetric_difference(set(R[l]))
                R[j] = sorted(list(symdif))
                if R[j] == []:
                    return (j)
                low_one[j] = R[j][0]       
    return (-1)

# Jeg har ingen store nye forbedringer for øyeblikket.
def general_dijkstra(cols,w,t,w_t=0):
    """
    cols -  An array of sets of indices representing the rows where the column has a one
    w    -  Array of weights of each column (len(w)=len(cols))
    t    -  Target vector: A set of indices
    w_t  -  Weight of the target (after adding columns that has to be in the solution)
    """
    
    (start,end) = (t.copy(),set())   # Starter i t, ender i 0.
    queue = []                       # Lager en tom liste for køen
    current = [start,w_t,[]]         # Initialbetingelser til algoritmen
    counter = 0                      # Antall gjennomkjøringer av while
    ub = float('inf')                # Upper bound. Kan kanskje regne ut et grovt estimat før.
    
    # Kjører så lenge den nåværende noden ikke er end-noden
    while (current[0]==end)==False:
        counter= counter +1
        # Nullstiller check som sier om vi har kolonner med vekt som går over upper bound fra current
        weight_check = False
        # Printer informasjon
        #if counter%100==0:
            #print("Current weight: ",current[1],"Step: ", counter, "Queue: ", len(queue), "Cols: ",len(cols))
        
        # Forced Path, neste rad som må ha en ener og kolonnene som har enere i raden
        (fp_row, next_cols) = face_choice(current,cols)   #fp_row = list(current[0])[0] er et naivt valg
        
        for i,c in next_cols:
            node = c.symmetric_difference(current[0])         # Den nye noden etter å summe på c på current
            weight = w[i] + current[1]                        # Den nye vekten
            
            # Ser om vekten kommer over upper bound, vi husker at det skjer og ser på neste simplex
            if weight > ub:
                weight_check = True
                continue
            
            path = current[2].copy()        # Stien opp til current
            path.append(c)                  # Stien til den nye noden
            
            # Check if there's a cycle in the path so far, if it is it is not a solution
            if cycle_in_path(path)!=-1:
                continue
            
            # Skjekker om noden allerede er i køen, og erstatter hvis den nye vekten er mindre enn den gamle.
            in_queue = False
            updating = False
            for j,u in enumerate(queue):
                if u[0] == node:
                    if weight < u[1]:
                        queue[j] = [node,weight,path]
                        updating = True
                    in_queue = True
                    break
            if in_queue == False:
                queue.append([node,weight,path])
                updating = True
            
            # Oppdaterer upper bound og køen hvis vi har oppdatert end-noden
            if updating == True:
                if node == end:
                    ub = weight
                    queue = sorted(queue,key=lambda x:x[1])       # Sorterer køen etter vekt
                    end_index = queue.index([node,weight,path])   # Finner index til end
                    queue = queue[:end_index+1]                   # Fjerner alt tyngre enn end
                    #print('UPPER BOUND UPDATED TO', ub)           # Printer den nye upper bound
        
        # Hvis vekten på en kolonne er for stor så kan vi fjerne den og alle større kolonner
        if weight_check == True:
            toobig = ub - current[1]   # Maksimalt hvor stor vekten kan være   
            for i in reversed(range(len(w))):
                if w[i] >= toobig:
                    cols.pop(i)     # Fjerner kolonner større enn maksimal
                    w.pop(i)        # Fjerner også den tilsvarende vekten
                    
        # Velger den neste (minste) node å se på
        if queue != []:
            minval = float('inf')
            for j,u in enumerate(queue):
                if u[1]<minval:
                    index = j
                    minval = u[1]  
            current = queue.pop(index)          # Velger den minste noden og fjerner den fra listen
    #print("Current weight: ",current[1],"Step: ", counter, "Queue: ", len(queue), "Cols =", len(cols))
    #print(t)
    return (current)   
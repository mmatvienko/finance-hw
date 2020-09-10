def generate_elementary_prices(srl, periods=10):
    # srl is the short rate lattice
    P = [[1]]

    def p_k_s(k, s):
        # find some elementary security at time k+1 and state s
        # going to the right anywhere but bottom
        first = P[k-1][s-1]/(2*(1 + srl[k-1][s-1]))
        second = P[k-1][s]/(2*(1 + srl[k-1][s]))
        return first + second
    
    def p_k_0(k):
        # going to the right at the bottom
        num = P[k - 1][0]
        den = 2*(1+srl[k - 1][0])
        return num / den
    
    def p_k_k(k):
        # going diagonally at the top
        num = P[k-1][k-1]
        den = 2*(1 + srl[k-1][k-1])
        return num / den
    
    for k in range(1, periods + 1):
        new_period = []
        
        # do bottom
        new_period.append(p_k_0(k))

        # do i - 1 middle items
        for s in range(1, k):
            new_period.append(p_k_s(k, s))
        
        # do top
        new_period.append(p_k_k(k))
        
        P.append(new_period)
        
    return P


def reduce_call(values, underlying, short_rate_lattice, n=6, K=80, qu=.5, qd=.5, american=True, debug=False):
    for period in range(n)[::-1]:
        # get the necessary information for this period
        curr_underlying = underlying[period]
        curr_short_rates = short_rate_lattice[period]
        previous_values = values[0]
        
        new_results = []
        
        for state in range(period + 1):
            discount_rate = 1/(1 + curr_short_rates[state])
            risk_neutral_val = (qu*previous_values[state] + qd*previous_values[state + 1])*discount_rate
            
            if american:
                val = max(max(curr_underlying[state] - K, 0), risk_neutral_val)
            else:
                val = max(0, risk_neutral_val)

            new_results.append(val)
        
        values.insert(0, new_results)
    return values

def reduce_future(future_lattice, periods, qu=.5, qd=.5):
    
    for period in range(periods)[::-1]:
        prev_prices = future_lattice[0]
        new_prices = []
        for state in range(len(prev_prices) - 1):
            new_price = qu*prev_prices[state] + qd*prev_prices[state + 1]
            new_prices.append(new_price)

        future_lattice.insert(0, new_prices)

    return future_lattice


def reduce_forward(forward_lattice, short_rate_lattice, periods, qu=.5, qd=.5, debug=True):
    
    # go through each time period
    for period in range(periods)[::-1]:
        new_prices = []
        prev_prices = forward_lattice[0]
        if debug: print(f"Period {period}\tlen curr price: {len(prev_prices) - 1}")
        if debug: print(f"Previous prices: {prev_prices}")

        # go through each state for the given time period
        for state in range(len(prev_prices) - 1):
            if debug: print(f"\t1st: {period}\t2nd: {state}")
                
            discount_rate = 1/(1 + short_rate_lattice[period][state])
            ex_discount_rnp = qu*prev_prices[state] + qd*prev_prices[state + 1]
            new_prices.append(discount_rate*ex_discount_rnp)

        forward_lattice.insert(0, new_prices)

    return forward_lattice

def generate_zcb_lattice(face_value, periods, short_rate_lattice, qu=0.5, qd=0.5, debug=False):
    # create values at maturity (face_value)
    lattice = [[face_value]*(periods + 1)]
    
    # step backward in time
    for i in range(periods):
        prev_prices = lattice[-1 - i]
        curr_short_rate = short_rate_lattice[periods - i  - 1]
        period_prices = []
        if debug: print(f"Period {periods - i  - 1}")
       
        # go through each state, there are x + 1 states in the xth period
        for state in range(periods - i):
            discount_rate = 1 / (1 + curr_short_rate[state])
            zij = discount_rate*(qu*prev_prices[state] + qd*prev_prices[state + 1])
            period_prices.append(zij)
            if debug: print(f"\tstate {state}\tshort rate: {curr_short_rate[state]}\tfirst pp: {prev_prices[state]}\tsecond pp: {prev_prices[state + 1]}")
        
        lattice.insert(0,period_prices)
        if debug: print(lattice) 

    return lattice


def generate_short_rate_lattice(r00, n, u, d):
    short_rate = [[r00]]

    for i in range(n):
        # get the latest
        latest_rates = short_rate[-1]
        new_rates = []

        # go through each short rate and scale it up
        for j in range(len(latest_rates)):
            new_rate = u * latest_rates[j]
            new_rates.append(round(new_rate,5))

        # only need one scale down since it "propogates"
        last = d * latest_rates[-1]
        new_rates.append(round(last,5))

        # add this periods rates to the lattice
        short_rate.append(new_rates)

    return short_rate


def reduce_swap(fixed_rate, short_rate_lattice, start=1, periods=10, qu=.5, qd=.5, debug=True):
    # once start time is achieve, just reduce using RNP
    # at the very end, only coupon is calc'ed 
    # in the middle, both coupon and underlying (RNP) are combined

    def coupon(period, state):
        # base on some period and state, gets the coupon payment
        if debug: print(f"period: {period} state: {state}")
        rate = short_rate_lattice[period][state]
        coupon = (rate - fixed_rate)/(1 + rate)

        return coupon
        
    if debug: print("started at the end")
    # build the final values that will be reduced
    swap_lattice = []
    for state in range(periods + 1):
        swap_lattice.append(coupon(periods, state))

    swap_lattice = [swap_lattice]
    if debug: print(swap_lattice[0])

    if debug: print("started middle")
    # reduce the swap intermediate part of the swap
    for period in range(start, periods)[::-1]:
        prev_coupon = swap_lattice[0] 
        new_coupon = []
    
        
        for state in range(period + 1):
            rate = short_rate_lattice[period][state]
            discount_rate = 1/(1 + rate)
            coupon = (rate - fixed_rate)
            
            rnp = (coupon + qu*prev_coupon[state] + qd*prev_coupon[state+1])*discount_rate
            new_coupon.append(rnp)
            
        swap_lattice.insert(0, new_coupon)
        if debug: print(swap_lattice[0])

   
    if debug: print("started the beginning part")
    # fill in the end if there's anything to do
    for i in range(start)[::-1]:
        prev_coupon = swap_lattice[0] 
        new_coupon = []

        for state in range(i + 1):
            discount_rate = 1/(1 + short_rate_lattice[i][state])
            rnp = (qu*prev_coupon[state] + qd*prev_coupon[state+1])*discount_rate

            new_coupon.append(rnp)

        swap_lattice.insert(0, new_coupon)
        if debug: print(swap_lattice[0])


    return swap_lattice

# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 17:04:13 2019

@author: Jared
Questions on research papers for Dr. Shindell
"""

'''
Tebaldi paper
Although we do not attempt to translate extreme heat degree days to yield impacts, we note that other authors have
assumed that yield loss increases proportional to accumulation of extreme heat measures when
plants are blooming (Deryng et al. 2014).
This might be difficult, but I wonder if it would be possible to accuont for extreme weather's impact on worldwide yields in some way?

The C4 crop, maize,
continues to show significant gains from mitigation, if of smaller size, while CO2 fertilization
of the C3 crop, wheat, generally cancels out the effects of mitigating climate changes in terms
of mean yield impacts (ignoring yield effects of exceeding 34 °C).

-- so this is saying that C4 crops are negatively effected by CO2, and C3 crops generally respond favorably to higher levels of CO2?
--I don't think it would be challenging to differentiate between these crops (because the only signifiant C4 crops is corn--beans are C3, right?) but I didn't see any actual coefficients we could use in the code for Corn
--I don't see where temperature threshold effects are addressed: they say that temperatures over 34 and 35 degrees celsius, at critical times in the year, can kill the reproductive processes in wheat and corn. I suppose we could look at the climatic variability of different regions around the world, and try to estimate the increase in # of days above 34 / 35 C (that's done in the paper), but I don't know how veracious that is. For example, what does "kill" mean? and what are those "windows in hte reproductive cycle", how long are they, when are they, do they vary according to growing season, how likely are they to stunt growth, how much yield loss does "stunt" entail, etc. 

'''

'''
MYERS ET. AL.
This says nothing about whether the relationship is linear: I.e we have no way of saying factually that an increase to say, 475 ppm, is proportionately less than 546. It just gives us the baseline and the change. 
Are we justified in assuming a linear response? Why?


Here we report that C3 grains and
legumes have lower concentrations of zinc and iron when grown
under field conditions at the elevated atmospheric CO2 concentration predicted for the middle of this century. C3 crops other than
legumes also have lower concentrations of protein, whereas C4
crops seem to be less affected.

We found that elevated [CO2] was associated with significant decreases in the concentrations of zinc and iron in all C3 grasses and legumes (Fig. 1 and Extended Data Table 1). For example, wheat grains
grown at elevated [CO2] had 9.3% lower zinc (95% confidence interval
(CI)212.7% to25.9%) and 5.1% lower iron (95% CI26.5% to23.7%)
than those grown at ambient [CO2].We alsofound that elevated [CO2]
was associated with lower protein content in C3 grasses, with a 6.3%
decrease (95% CI 27.5% to 25.2%) in wheat grains and a 7.8% decrease
(95% CI 28.9% to 26.8%) in rice grains. Elevated [CO2] was associated
with a small decrease in protein infield peas, and therewas no significant
effect in soybeans or C4 crops (Fig. 1 and Extended Data Table 1).

'''

'''
ZHOU ET. AL

QUESTION: the numbers they got for temperature changes are, in some cases, very different from our coefficients. What gives? 

QUESTION: I think this conclusion is missing half of the fact. It seems pretty obvious that
wheat yields at 40 degrees farenheit won't be higher than at 80 for the simple reason that colder sites have shorter growing seasons and thus lower yields. There must be some "sweet spot" in the total yield
and this sweet spot, I assume, must not only be impacted by the limitations of growing seasons but also by the fact (this is a big assumption, and it is half of my question) that biological processes within crops probably slow at lower temperatures. Or at least there's a point at which the negative effects of temperature rise stop.
Regardless, it's a fact that colder climate agricultural regions like the upper midwest, southern canada, and russia should expect to see, and in some cases already have seen, large increases in crop yield. 
I think we'd be remiss if we didn't address this as well. 


 Without CO2 fertilization,
effective adaptation, and genetic improvement, each degree-Celsius
increase in global mean temperature would, on average, reduce
global yields of wheat by 6.0%, rice by 3.2%, maize by 7.4%, and
soybean by 3.1%. Results are highly heterogeneous across crops
and geographical areas, with some positive impact estimates. Multimethod analyses improved the confidence in assessments of future
climate impacts on global major crops and suggest crop- and regionspecific adaptation strategies to ensure food security for an increasing world population.
'''

'''
Deutsch.

QUESTION: I'm skeptical about this one right off the bat: Pick a corn farm in Iowa and go 500 miles south and pick another corn farm in Arkansas or Texas. Will this second farm lose a higher proportion of its yield to insects than the first? 
It's warmer there, and we shouldn't forget that the very simplest way to model how an area will change due to global warming, in the long run, is to move south, or down to a lower elevation. 
I hope the paper considers this. 

I think it seems to, it suggests that there are other ecological factors at play here, but that the dominant trend of "more insects, higher temperatures, more plants being eaten by insects" holds quite broadly. 

QUESTION: I notcied one of the figures in here features a crop model that has been clipped so that only regions that produce the crop in question are colored in. I need to figure out how to do that. Do you know?
Also I think it's a little difficult to deal with / dubious that their units are in change in yield loss rather than % loss as a whole. Feels like they're trying to conceal how small their numbers actually are. 


When average global surface temperatures
increase by 2°C, the median increase in yield
losses owing to pest pressure is 46, 19, and 31%
for wheat, rice, and maize, respectively, bringing
total estimated losses to 59, 92, and 62 metric
megatons per year


'''
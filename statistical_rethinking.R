library(rethinking)


#Ch. 2 Important Code

#Grid Approximation Code

#something that I think I've basically not thought through properly is the fact
#that what we are doing with this kind of estimation is giving probabilities
#to probability values e.g what happens here is that we are trying to 
#estimate the true probability of landing on water - and we do this by assigning
#probabilities to all of the different p values. Roughly, I think what happens
#below is that "p_grid" is basically a grid containing 20 possible values for
#the probability of water. Then, a prior probability is assigned to each of these
#values i.e it is decided how likely each value will be. Then we calculate the
#likelihood by using the binomial distribution to determine how likely it is
#that each p value produces the data, where the relevant data is implicit in
#the binomial calculation i.e. the data is that there was 9 trials and on 6 
#water was found. Importantly for your intuitions, the likelihood calculation
#doesn't involve the prior at all (I think I was confused and did think that
#something along those lines was occurring, but apparently not). All the likelihood
#is is the calculation of how likely it is that the actual data we have would occur
#for a given P value. How likely it is that the given P value actually obtains
#is what the prior and posterior probabilities are concerned with. Binomial is used
#for this because there are only two possible outcomes for each interation. In the
#code, what happens is that each p value is pumped in, and the binom outputs
#the likelihood that the data would emerge on that p value, for each p value in
#the prior distribution.
#
#We then need to actually work on calculating the posterior probability. We do
#this by multiplying the prior probability (which indexes the likelihood that
#any given P value actually obtains) by the likelihood (which indexes the likelihood
#that the data we have obtained would emerge from any of the p values). I'm trying to get
#a proper intuition on this -- for an individual p, we have the prior likelihood that this p
#is true, and the likelihood that the data we have would emerge if this p was true.
#Imagine a case where our prior was that P1 was a really likely p value but that the 
#likelihood value of P1 showed that the data was really unlikely if P1 obtained. 
#Multiplying them will thus alter the probability of P1 by conditioning on the data 
#(where the data indicates that P1 is unlikely), and therefore decreasing the probability
#of P1 in the posterior distribution. In this way, the likelihood, which represents the data
#essentially disciplines the prior distribution by increasing the probability of p values
#which are likely to create the data/evidence and decreasing the likelihood of p values which are 
#unlikely to create the data evidence. This is the essence of the bayesian update.
#
#Then we need to understand the purpose of the regularising variable or marginal/average
#likelihood. On the surface this is just needed to make sure that all of the probabilities
#in the posterior distribution sum to one. What the value actually is is the average of 
#the likelihoods of producing the data/evidence for each p being considered.
#It might be more intuitive to think of this in terms of the bayes equation being 
#prior probability * likelihood//marginal likelihood. This makes it clear that what happens
#is that the likelihood is regularised into a proportion of the total likelihood across all
#p values. What the marginal likelihood reflects more precisely is the sum of all priors
#with their respective likelihoods. Why this is referred to as the "average" or "marginal" likelihood
#I am not so sure. The most important thing seems to be that when you divide by this you - obvs -
#get the result that all of the individual likelihood + prior combinations collectively sum to 1. 
#
#Now to explain the actual grid approximation we are doing here. Basically we prespecify our 20
#possible p values instead of doing an infinite calculation across all of the p values. It's clear
#that doing the infinite calculation would soon make calculating the marginal likelihood 
#(and likelihood? not sure) intractible. With grid approximation, all we do is sum across our
#likelihood/probability combos to create the marginal likelihood, and then divide each individual
#combo by this value in order to regularise the posterior distribution and ensure that it's 
#probabilities sum to one. 
#Apparently this is a relatively suboptimal way to go
#about things but it makes the intuitions behind what we are doing clear. 



# define grid
p_grid <- seq( from=0 , to=1 , length.out=20 )
# define prior
#prior <- exp( -5*abs( p_grid - 0.5 ) )
# and some alternative priors to try:
#prior <- ifelse( p_grid < 0.5 , 0 , 1 )
#prior <- exp( -5*abs( p_grid - 0.5 ) )

# compute likelihood at each value in grid
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )


#Quadratic Approximation

#I'm not sure on the details of this - they aren't really explained - but the
#essence is that grid approximation because extremely costly when there are lots
#of parameters (because calculations increase exponentially with more parameters)
#Quadratic Approx - which basically just assumes that the posterior distribution
#will be gaussian - gets around this because the distribution an be completely
#described using only a mean and variance. These are calculated using an optimisation
#algorithm to locate the mode (highest probability p value) of the posterior distribution,
#and then estimate the curvature around the peak to derive the quadratic for the whole
#distribution. Here is an example -- it's all explained more later apparently. 
#
#It uses the "map" tool from the rethinking package. You must input a list of 
#starting values for the likelihood and prior parameters, as well as a list of data.
#It outputs a mean and standard deviation which defines a gaussian posterior distribution.



library(rethinking)
globe.qa <- map(
  alist(
    w ~ dbinom(9,p) , # binomial likelihood
    p ~ dunif(0,1) # uniform prior
  ) ,
  data=list(w=6) )
# display summary of quadratic approximation
precis( globe.qa )

#Apparently the following "analytical calculation" produces the exact posterior
#without estimation. This is not explained but can apparently be easily looked up.
# The code below plots both the quad approx and analytic curves so you can see
#the weaknesses of quad approach. The quad approx will get increasingly better the
#more data you pump in (better here meaning closer to the exact posterior calculated
#by the analytical approach, better does NOT mean more accurate necessarily
#just more optimal given assumptions. Ultimately data determines the horizon of
#possible accuracy). We are still here just trying to approximate the ideal calculations
#for the data that we have. This also emphasises that much of the time when stat
#procedures want more data, it's often not to maximise real world accuracy but
#rather estimation-of-bayes accuracy, which is clearly distinct through also 
#presumably related. Interesting!


# analytical calculation
w <- 6
n <- 9
curve( dbeta( x , w+1 , n-w+1 ) , from=0 , to=1 )
# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )

#MCMC

#basically there are some important model types like multi-level mixed effects models
#for which grid approx or quadratic approx are not always satisfactory. Insane
# No of params + mysterious reasons why quad doesn't work either in some cases. 
#MCMC brought in in such cases. Instead of attempting to compute the post dist,
#it merely draws samples from the post dist. You then build a picture of the post
#dist using these samples. Normally we work directly with the samples instead of 
#estimating the actual post dist from them. Working with samples is focus of Next Ch.

#Ch.3 Sampling the Imaginary 

#Classic bayes e.g. -- vamp test, detects vampirism 95% of the time - 
#Pr(positive|vamp) = 0.95. False positive 1% of the time so Pr(positive|mortal) = 0.01.
#And vamps are only 0.1% of the total pop so Pr(vamp) = 0.001. What is prob of vamp
#given positive test? Can use bayes. Pr(positive) is just the average probability of
#a positive test result. Pr(pos) = Pr(pos|vamp)Pr(vamp) + Pr(post|mort)(1-Pr(vamp)).


PrPV <- 0.95
PrPM <- 0.01
PrV <- 0.001
PrP <- PrPV*PrV + PrPM*(1-PrV)
( PrVP <- PrPV*PrV / PrP ) # 0.0868

#Whenever a condition of interest is very rare, a test that finds all true cases
#is useless if it's relative false positive rate is high enough, as it results in
#a signal far too noisy to act upon or to provide any meaningful information.

#Thinks this e.g. is bad because it doesn't exhibit distinctly bayesian qualities.
#Probabilities here are assigned to events, not theoretical hypotheses which drive events
#and thus no one disagrees with bayes theorem here. It also makes bayes seem hard. 
#Notes that if you present this problem in terms of counts (natural frequencies)
#it becomes much more intuitive. Whatever the cause, this is now exploited by sampling
#from distributions. Post dist is still a prob dist that we an plausibly draw samples
#from. The sampled events are parameter values (of the distribution we are sampling)
#The post dist defines the expected frequency that different parameter values will appear
#once we start plucking parameters out of it

#Sampling approach is good because it lets you dodge integrals whilst still getting 
#the same info. Sampling from the posterior ("an empirical attack") lets you ask Qs
#of a model without the advanced mathematics. And MCMC etc produce nothing but samples anyways.
#Fun note showing that if you assume that the general probability of uncovering the
#truth is very small, then for 0.95 testing the likelihood that the hypothesis under test
#is actually true is unlikely. Shows via bayes that we must be very careful even with
#our best testing mechanisms. 

#generating some samples -- redo globe tossing grid approx

p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

#Now we will sample from the posterior. This is just like sampling from IRL - 
#the post dist values near the peak are more common than those near the tails (gaussian)
#just like how things more common in the world will be more commonly sampled. 
#If one samples extensively from the distribution, the distribution of samples should
#just map quite well onto the distribution of the world (again, like normal sampling)

#Getting 10000 samples from the distribution. sample() randomly pulls values from a vector
#Vector here is p_grid, probability is from the post dist. In essence, this is just
#estimating a distribution we already created with grid approx

samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
plot(samples)
dens(samples)

#Common Qs for a post dist: 
#         (1) How much posterior probability lies below some parameter value?
#         (2) How much posterior probability lies between two parameter values?
#         (3) Which parameter value marks the lower 5% of the posterior probability?
#         (4) Which range of parameter values contains 90% of the posterior probability?
#         (5) Which parameter value has highest posterior probability?
#

#Qs about intervals of defined boundaries

#probability that the proportion of water is < 0.5. for grid approx post dist
#you can just add up all the probabilities where the corresponding parameter value
#is < 0.5

# add up posterior probability where p < 0.5
sum( posterior[ p_grid < 0.5 ] ) # 0.1719

#I've lost my intuition about what is happening. The post dist represents different
#values for how much of the Earth is covered in water. The various percentages are
#assigned probabilities in the actual post dist. When we sample the post dist a lot
#we basically copy it. But for now what we are dealing with is the actual posterior. 
#p_grid is the grid of water probabilities. So what we are doing
#is summing the probabilities assigned to water percentages below 0.5 (sums to 0.1719)

#we now want to do the same with the samples we just got, because that's more generalisable.
#each element of the samples is a percentage for how much of the earth is covered in water.
#unlike in the actual posterior, the probabilities of each percentage (from p_grid)
#is not directly represented. So what we need to do is literally count the amount of
#values in the sample that are below 0.5 and divide this by the total number of values
#in order to get a percentage. This is here done again using the "sum" function because 
#of a quirk with how sum() works, but this can seem unintuitive. In cases like this, sum()
#functions like table()--where the latter is a straightforward counting function-- and this
#is probably a less confusing way to represent what happens. 

sum( samples < 0.5 ) / 1e4
#and you get the same output with table()
table(samples<0.5)

#we can also do this to find out about other regions of the post dist

sum( samples > 0.5 & samples < 0.75 ) / 1e4

#Qs relating to intervals of defined mass

#confidence intervals basically. also or should be called "credible interval" when working with
#post dist.They report two parameter values that contain between them a specified amount
#of post probability - probability mass. Easier to do this with sample rather than actual dist.
# e.g. want to find the lower 80% probability. You know it has to start at 0
quantile( samples , 0.8 )
#or the middle 80% - between 10th and 90th percentile
quantile( samples , c( 0.1 , 0.9 ) )

#these percentile intervals can do a good job of communicating the shape of the
#distribution, as long as the distribution is not too asymmetrical. Basically, 
#these communicate information under assumptions of a gaussian. But if it's not
#gaussian, chopping out the middle might be very misleading. The middle might not
#be the densest or most important part of the dataset. e.g. in this posterior --
#we get the percentile intervals, with 25% of probability mass above and below
#the interval. But this hides that the top 25% is taken up by a small amount of 
#highly probably values.

p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep(1,1000)
likelihood <- dbinom( 3 , size=3 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid , size=1e4 , replace=TRUE , prob=posterior )

PI( samples , prob=0.5 )

#we can find the narrowest interval containing the specified probability mass 
#(in other words, the relevant peak) by finding the highest posterior density interval
HPDI( samples , prob=0.5 )

#mostly these two intervals looks the same if the post dist isn't v skewed.
#HDPI is more computationally intensive and has greater simulation variance i.e. is
#sensitive to how many samples are drawn from the post dist. but mostly it shouldn't
#matter that much. 









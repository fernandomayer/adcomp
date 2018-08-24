library(TMB)
compile("thetalog2.cpp")
dyn.load(dynlib("thetalog2"))

## Read data
Y <- scan("thetalog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y)
plot(Y)

## Parameter initial guess
parameters <- list(
  X = data$Y*0,
  logr0 = 0,
  logtheta = 0,
  logK = 6,
  logQ = 0,
  logR = 0
)

## Fit model
## Inicialmente vou fazer sem assumir que X é aleatório, para ver se a
## verossimlhança será a mesma.
obj <- MakeADFun(data, parameters,
                 ## random="X",
                 DLL="thetalog2")
obj$fn()
## Como chegar nesse valor?

## Vou fazer uma função que tenta "imitar" fielmente o programa escrito
## em C++ (adaptando algumas coisas devido à diferença de linguagens).
## Veja no artigo do TMB no JSS (pg. 13) que a função de verossimilhança
## conjunta é a soma das duas individuais (de X e de Y).
thetalog <- function(Y, X, logr0, logtheta, logK, logQ, logR){
    r0 <- exp(logr0)
    theta <- exp(logtheta)
    K <- exp(logK)
    Q <- exp(logQ)
    R <- exp(logR)
    timeSteps <- length(Y)
    m <- numeric(length(timeSteps))
    ans <- numeric(length(timeSteps))
    for(i in 2:timeSteps){
        ## Note que o X é utilizado para o cálculo da verossimilhança,
        ## mas nesse caso, X[i-1] e X[i] serão sempre 0
        m[i-1] <- X[i-1] + r0*(1 - (exp(X[i-1])/K)^theta)
        ans[i-1] <- dnorm(X[i], m[i-1], sqrt(Q), log = TRUE)
    }
    ans <- -sum(ans) # vero de X
    ans2 <- numeric(length(timeSteps))
    for(i in 1:timeSteps){
        ## Note que o X[i] é sempre 0
        ans2[i] <- dnorm(Y[i], X[i], sqrt(R), log = TRUE)
    }
    ans2 <- -sum(ans2) # vero de Y
    res <- ans + ans2  # vero conjunta
    ## Para "recompor" o X. Note que se eu fizesse isso antes, iria
    ## alterar a verossimilhança do Y.
    for(i in 2:timeSteps){
        X[i] <- X[i-1] + r0*(1 - (exp(X[i-1])/K)^theta)
    }
    saida <- list(m = m, X = X, nllX = ans, nllY = ans2, nll = res)
    return(saida)
}

## Aplica a função
m0 <- thetalog(Y = Y, X = numeric(length(Y)),
               logr0 = 0, logtheta = 0, logK = 6,
               logQ = 0, logR = 0)
m0[3:5]

## Portanto note que são iguais
m0$nll
obj$fn(obj$par)

## Para ver o que é o X (estados)
range(Y)
ylim <- c(0, 8)
plot(Y, ylim = ylim)
lines(m0$X)

## E veja que o m é igual pra todo mundo
sum(diff(m0$m))

## É claro que do jeito que está, essa função não pode ser otimizado
## pois X está sendo tratado como paraâmetro, então existe um parâmetro
## para cada observação.

## Por isso, quando usamos o argumento random = "X" na MakeADFun,
## estamos dizendo que o X são os efeitos aleatórios, e que devem ser
## "integrados fora" da função de verossimilhança (usando a aproximação
## de Laplace). Veja o que diz em ?MakeADFun:

## Random effects are specified via the argument 'random': A component
## of the parameter list is marked as random if its name is matched by
## any of the characters of the vector 'random' (Regular expression
## match is performed if 'regexp=TRUE'). If some parameters are
## specified as random effects, these will be integrated out of the
## objective function via the Laplace approximation. In this situation
## the functions 'fn' and 'gr' automatically perform an optimization of
## random effects for each function evaluation. This is referred to as
## the 'inner optimization'. Strategies for choosing initial values of
## the inner optimization can be controlled via the argument
## 'random.start'. The default is 'expression(last.par.best[random])'
## where 'last.par.best' is an internal full parameter vector
## corresponding to the currently best likelihood. An alternative choice
## could be 'expression(last.par[random])' i.e. the random effect
## optimum of the most recent - not necessarily best - likelihood
## evaluation. Further control of the inner optimization can be obtained
## by the argument 'inner.control' which is a list of control parameters
## for the inner optimizer 'newton'.

## Assim, precisamos declarar X como random, e assim podemos prosseguir
## com a otimização
obj2 <- MakeADFun(data, parameters,
                  random="X",
                  DLL="thetalog2")
obj2$fn(obj2$par)
## Note que é feita uma otimização para se chegar nesse valor. Esse é o
## método de Newton (inner optimization) usado para fazer a aproximação
## de Laplace (calcular a integral na dimensão de X). O resultado é
## claramente diferente de
obj$fn(obj$par)

## Agora podemos prosseguir com a otimização
system.time(opt <- nlminb(obj2$par, obj2$fn, obj2$gr))
rep <- sdreport(obj2)
rep
rep$par.fixed
rep$par.random

## Para ver o que é o X (estados)
range(Y)
ylim <- c(0, 8)
plot(Y, ylim = ylim)
lines(m0$X)
lines(rep$par.random, col = 2)
## \o/

## VER:

## - Como fazer uma predicao
## - Mostrar que a funcao thetalog pode usar os argumentos sem serem
## passados em log
## - Tentar escrever uma funcao (parecida com thetalog) que possa ser
## otimizada direto pela optim. Talvez assumindo erro somente nas
## observacoes para ficar so com uma verossimilhança. E tentar fazer uma
## verossimilhança tambem no TMB

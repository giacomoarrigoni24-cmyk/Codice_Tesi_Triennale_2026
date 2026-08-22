# La funzione legame è un oggetto che permette di organizzare i diversi legami implementati per i dati binari. Risulta importante precisare come sia possibile aggiungere altre funzioni legame ed utilizzarle poi per le altre funzioni presenti nella repository senza dover effettuare altre modifiche.

.links_GLM <- function() {
  list(
    logit = list(
      G     = function(x) exp(x) / (1 + exp(x)),
      g     = function(x) (exp(x) / (1 + exp(x))) * (1 / (1 + exp(x))),
      gdash = function(x) (exp(x) * (1 - exp(x))) / (1 + exp(x))^3,
      invG  = function(x) log(x / (1 - x))
    ),
    probit = list(
      G     = pnorm,
      g     = dnorm,
      gdash = function(x) -x * dnorm(x),
      invG  = qnorm
    ),
    cloglog = list(
      G     = function(x) 1 - exp(-exp(x)),
      g     = function(x) exp(x - exp(x)),
      gdash = function(x) (1 - exp(x)) * exp(x - exp(x)),
      invG  = function(x) log(-log(1 - x))
    ),
    loglog = list(
      G     = function(x) exp(-exp(-x)),
      g     = function(x) exp(-x - exp(-x)),
      gdash = function(x) (-1 + exp(-x)) * exp(-x - exp(-x)),
      invG  = function(x) -log(-log(x))
    ),
    cauchit = list(
      G     = function(x) 0.5 + atan(x) / pi,
      g     = function(x) (1 / pi) / (1 + x^2),
      gdash = function(x) -(2 * x / (1 + x^2)) * (1 / pi) / (1 + x^2),
      invG  = function(x) tan((x - 0.5) * pi)
    )
  )
}


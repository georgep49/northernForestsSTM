logistic <- function(r, K, x) {K / (1 + exp(-r*x))}

plot(logistic(r = 0.01, K = 1.0, x = seq(0,1,0.01)))


# https://math.stackexchange.com/questions/4620969/logistic-function-from-0-to-1-for-input-0-to-1
fx <- function(x, x0, k) {1 / (1 + ((1-x) * x0 / x * (1 - x0)) ^ k)}

plot(fx(x = seq(0, 1, 0.01), x0 = 0.5, k = 2))


## declining and increasing funcs

incr_limit <- function(b1, b2, x) {b1 * (1 - exp(-b2 * x))}
decr_limit <- function(b1, b2, b0, x) {b1 *  exp(-b2 * x) + b0}

y <- incr_limit(b1 = 1, b2 = 0.1, x = seq(0, 500, 5) / 10)
y <- decr_limit(b1 = 1, b2 = 0.1, b0 = 0, x = seq(0, 500, 5) / 10)

1 - ((1 / 1:8) ^ 0.2)

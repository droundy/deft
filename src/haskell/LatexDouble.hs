{-# LANGUAGE PatternGuards #-}

module LatexDouble where

latexDouble :: Double -> String
latexDouble 0 = "0"
latexDouble x | x < 0 = "-" ++ latexDouble (-x)
latexDouble x
  | Just (n,d,p) <- double2powfrac x 1 = latexPowRat (n,"",d,"",p)
  | Just (n,d,p) <- double2powfrac x pi = latexPowRat (n,"\\pi",d,"",p)
  | Just (n,d,p) <- double2powfrac x (1/pi) = latexPowRat (n,"",d,"\\pi",p)
  | Just (n,d,p) <- double2powfrac x (pi**1.5) = latexPowRat (n,"\\pi^{3/2}",d,"",p)
  | otherwise = show x

latexPowRat :: (Int, String, Int, String, Double) -> String
latexPowRat (n,ns,d,ds,1) = latexRat (n,ns,d,ds)
latexPowRat (n,ns,d,ds,0.5)
  | otherwise = "\\sqrt{" ++ latexRat (n,ns,d,ds) ++ "}"
latexPowRat (n,ns,d,ds,p) = "\\left(" ++ latexRat (n,ns,d,ds) ++ "\\right)^{" ++ latexDouble p ++ "}"

latexRat :: (Int, String, Int, String) -> String
latexRat (n,ns,1,"") = latexProd (n,ns)
latexRat (n,"",d,"")
  | n < 10 && d < 10 = "\\tfrac" ++ show n ++ show d
latexRat (n,ns,d,ds) = "\\frac{" ++ latexProd (n,ns) ++ "}{" ++ latexProd (d,ds) ++ "}"

latexProd :: (Int, String) -> String
latexProd (n,"") = show n
latexProd (-1,s) = "-" ++ s
latexProd (1,s) = s
latexProd (n,s) = show n ++ s

double2powfrac :: Double -> Double -> Maybe (Int,Int,Double)
double2powfrac x times = powby ([1..5] ++ [0.5,1.5 .. 4.5])
  where powby [] = Nothing
        powby (p:ps) = case double2frac (x ** (1/p) / times) of
                         Just (n,d) -> Just (n,d,p)
                         Nothing -> powby ps

double2frac :: Double -> Maybe (Int, Int)
double2frac f = do denominator <- divby [1..36]
                   numerator <- double2int (f * fromIntegral denominator)
                   Just (numerator, denominator)
  where divby (d:ds) = case double2int (f * fromIntegral d) of
                         Just _ -> Just d
                         Nothing -> divby ds
        divby [] = Nothing

double2int :: Double -> Maybe Int
double2int f = if abs(fromInteger (round f) - f) < 1e-13*f
               then Just (round f :: Int)
               else Nothing

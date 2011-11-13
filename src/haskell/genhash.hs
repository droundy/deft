import Control.Monad ( mapM_ )
import System.Random ( randomRIO )

main = do let allchrs = [' '..'~']
              goodchrs = ['a'..'z'] ++ ['A'..'Z'] ++ ['0'..'9']
              good = do n <- randomRIO (0, length goodchrs - 1)
                        return $ goodchrs !! n
              mix n = do putStr "mix "
                         putStr $ show n
                         putStr " = "
                         good >>= putStrLn . show
          putStrLn "\\begin{code}"
          putStrLn "module Hash ( hash ) where\n"
          putStrLn "import Data.Char ( ord )\n"
          putStrLn "hash :: String -> String"
          putStrLn "hash s = [hashc 'a' s, hashc 'b' s, hashc 'c' s, hashc 'd' s]"
          putStrLn "  where hashc c \"\" = c"
          putStrLn "        hashc c (x:xs) = hashc (combine c x) xs"
          putStrLn "        combine c x = mix (ord c + ord x)"
          putStrLn "\nmix :: Int -> Char"
          let nmax = show (length allchrs)
          putStrLn $ "mix n | n > " ++ nmax ++ " = mix (n `mod` " ++ nmax ++ ")"
          mapM_ mix [0..length allchrs]
          putStrLn "mix _ = 'x'"
          putStrLn "\\end{code}"

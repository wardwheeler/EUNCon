{- |
Module      :  GeneralUtilities.hs
Description :  Module with useful functions
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

module GeneralUtilities where

import           Data.Array
import qualified Data.Text  as T
import           System.Random
import           Data.Array.IO
import           Control.Monad
import           Data.Hashable
import           Data.List
import           Data.Time
import           Data.Time.Clock.POSIX
import           System.IO.Unsafe
import           Text.Read
import           Data.Maybe
-- import           Data.Global -- won't compile perhaps ghc-9 issue




-- | functions for triples, quadruples
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

snd3 :: (a,b,c) -> b
snd3 (_,e,_) = e

thd3 :: (a,b,c) -> c
thd3 (_,_,e) = e

fst4 :: (a,b,c,d) -> a
fst4 (e,_,_,_) = e

snd4 :: (a,b,c,d) -> b
snd4 (_,e,_,_) = e

thd4 :: (a,b,c,d) -> c
thd4 (_,_,e,_) = e

fth4 :: (a,b,c,d) -> d
fth4 (_,_,_,f) = f

fst5 :: (a,b,c,d,e) -> a
fst5 (e,_,_,_,_) = e

snd5 :: (a,b,c,d,e) -> b
snd5 (_,e,_,_,_) = e

thd5 :: (a,b,c,d,e)-> c
thd5 (_,_,e,_,_) = e

fth5 :: (a,b,c,d,e) -> d
fth5 (_,_,_,e,_) = e

fft5 :: (a,b,c,d,e) -> e
fft5 (_,_,_,_,e) = e

fst6 :: (a,b,c,d,e,f) -> a
fst6 (e,_,_,_,_,_) = e

snd6 :: (a,b,c,d,e,f) -> b
snd6 (_,e,_,_,_,_) = e

thd6 :: (a,b,c,d,e,f)-> c
thd6 (_,_,e,_,_,_) = e

fth6 :: (a,b,c,d,e,f) -> d
fth6 (_,_,_,e,_,_) = e

fft6 :: (a,b,c,d,e,f) -> e
fft6 (_,_,_,_,e,_) = e

six6 :: (a,b,c,d,e,f) -> f
six6 (_,_,_,_,_,e) = e

-- | editDistance is a naive edit distance between two lists
-- takes two  lists and returns edit distance
--- from  https://wiki.haskell.org/Edit_distance
editDistance :: Eq a => [a] -> [a] -> Int
editDistance xs ys = table ! (m,n)
    where
    (m,n) = (length xs, length ys)
    x     = array (1,m) (zip [1..] xs)
    y     = array (1,n) (zip [1..] ys)

    table :: Array (Int,Int) Int
    table = array bnds [(ij, dist ij) | ij <- range bnds]
    bnds  = ((0,0),(m,n))

    dist (0,j) = j
    dist (i,0) = i
    dist (i,j) = minimum [table ! (i-1,j) + 1, table ! (i,j-1) + 1,
        if x ! i == y ! j then table ! (i-1,j-1) else 1 + table ! (i-1,j-1)]

-- | getBestMatch compares input to allowable commands and checks if in list and if not outputs
-- closest match
-- call with (maxBound :: Int ,"no suggestion") commandList inString
getBestMatch :: (Int, String) -> [String] -> String -> (Int, String)
getBestMatch curBest@(minDist, _) allowedStrings inString =
    if null allowedStrings then curBest
    else
        let candidate =  head allowedStrings
            candidateEditCost = editDistance candidate inString
        in
        if candidateEditCost == 0 then (0, candidate)
        else if candidateEditCost < minDist then getBestMatch (candidateEditCost, candidate) (tail allowedStrings) inString
        else getBestMatch curBest (tail allowedStrings) inString

-- | getCommandErrorString takes list of non zero edits to allowed commands and reurns meaningful error string
getCommandErrorString :: [(Int, String, String)] -> String
getCommandErrorString noMatchList =
    if null noMatchList then ""
    else
        let (_, firstCommand, firstMatch) = head noMatchList
            firstError = "\tBy \'" ++ firstCommand ++ "\' did you mean \'" ++ firstMatch ++ "\'?\n"
        in
        firstError ++ getCommandErrorString (tail noMatchList)

-- | isSequentialSubsequence takes two lists and determines if the first List is
-- a subsequence of the second but the elements must be sequencetial unlike
-- isSubsequenceOf in Data.List
-- Uses Text.filter to see if there is a match
--isSequentialSubsequence :: (Eq a) => [a] -> [a] -> Bool
isSequentialSubsequence :: String -> String -> Bool
isSequentialSubsequence firstL secondL
  | null firstL = False
  | length firstL > length secondL = False
  | otherwise =
    let foundNumber = T.count  (T.pack firstL) (T.pack secondL)
    in
    foundNumber /= 0

-- | shuffle Randomly shuffles a list
--   /O(N)/
-- from https://wiki.haskell.org/Random_shuffle
shuffle :: [a] -> IO [a]
shuffle xs = do
        ar <- newArrayLocal  n xs
        forM [1..n] $ \i -> do
            j <- randomRIO (i,n)
            vi <- readArray ar i
            vj <- readArray ar j
            writeArray ar j vi
            return vj
  where
    n = length xs
    newArrayLocal :: Int -> [a] -> IO (IOArray Int a)
    newArrayLocal  nL xsL =  newListArray (1,nL) xsL



-- | randomList generates a random list from a seed--no IO or ST monad 
-- but needs a good seed
randomList :: Int -> [Double]
randomList seed = randoms (mkStdGen seed) :: [Double]


-- | selectListCostPairs is generala to list of (a, Double)
-- but here used for graph sorting and selecting)takes a pair of graph representation (such as String or fgl graph), and
--- a Double cost and returns the whole of number of 'best', 'unique' or  'random' cost
-- need an Eq function such as '==' for Strings or equal for fgl
-- optionsList must all be lower case to avoicd
-- assumes options are all lower case
-- options are pairs of Sting and number for number or graphs to keeep, if number is set to (-1) then all are kept
-- if the numToKeep to return graphs is lower than number of graphs, the "best" number are returned
-- except for random.
selectListCostPairs :: (Eq a, Hashable a) => (a -> a -> Bool) -> [(a, Double)] -> [String] -> Int -> [(a, Double)] 
selectListCostPairs compFun pairList optionList numToKeep = 
  if null optionList then error "No options specified for selectGraphCostPairs"
  else if null pairList then []
  else 
    let firstPass = if ("unique" `elem` optionList) then nubBy compFunPair pairList
                    else pairList
        secondPass = if ("best" `elem` optionList) then reverse $ sortOn snd firstPass
                     else firstPass
    in
    if ("random" `notElem` optionList) then take numToKeep secondPass 
    else -- shuffling with hash of structure as seed (not the best but simple for here)
      let seed = getSystemTimeSecondsUnsafe -- hash pairList
          randList = randomList seed
          pairListWRand = zip randList secondPass
          thirdPass = fmap snd $ sortOn fst pairListWRand

      in
      take numToKeep thirdPass
  where compFunPair fGraph sGraph = compFun (fst fGraph) (fst sGraph)

-- | getSystemTimeSeconds gets teh syste time and returns IO Int
getSystemTimeSeconds :: IO Int
getSystemTimeSeconds = do
    systemTime <- getCurrentTime  
    let timeD = (round $ utcTimeToPOSIXSeconds systemTime) :: Int
    return timeD

{-# NOINLINE getSystemTimeSecondsUnsafe #-}
-- | getSystemTimeSecondsUnsafe gets the system time and returns Int via unsafePerformIO
-- without the NOINLINE the function would probbaly be comverted to a 
-- constant which would be "safe" and OK as a random seed or if only called once
getSystemTimeSecondsUnsafe :: Int
getSystemTimeSecondsUnsafe = unsafePerformIO getSystemTimeSeconds
    
-- | stringToInt converts a String to an Int
stringToInt :: String -> String -> Int
stringToInt fileName inStr = 
  let result = readMaybe inStr :: Maybe Int
  in
  if result == Nothing then errorWithoutStackTrace ("\n\n'Read' 'tcm' format error non-Integer value " ++ inStr ++ " in " ++ fileName)
  else fromJust result

-- | makeIndexPairs takes n and creates upper triangular matrix pairs (0,m)
makeIndexPairs :: Bool -> Int -> Int -> Int -> Int -> [(Int, Int)]
makeIndexPairs doDiagValues numI numJ indexI indexJ =
    if indexI == numI then []
    else if indexJ == numJ then makeIndexPairs doDiagValues numI numJ (indexI + 1) 0
    else 
        if doDiagValues && (indexI == indexJ) then (indexI, indexJ) : makeIndexPairs doDiagValues numI numJ indexI (indexJ + 1)
        else if (indexI < indexJ) then (indexI, indexJ) : makeIndexPairs doDiagValues numI numJ indexI (indexJ + 1)
        else makeIndexPairs doDiagValues numI numJ indexI (indexJ + 1)
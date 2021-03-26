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
import qualified Data.Text       as T


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








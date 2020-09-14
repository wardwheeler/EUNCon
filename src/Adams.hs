{- |
Module      :  adams.hs 
Description :  functions to create Adams II consensus trees
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-# LANGUAGE DeriveAnyClass, DeriveGeneric #-}

module Adams (makeAdamsII) where 

import Data.List
import Data.Maybe
import qualified Data.Vector as V
import qualified Data.Graph.Inductive.Graph as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import Debug.Trace
import qualified Data.Set as Set
import qualified Data.BitVector as BV
import Control.Parallel.Strategies
import Control.Concurrent
import System.IO.Unsafe


{-
-- NFData instance for parmap/rdeepseq
instance NFData BV.BV where
  rnf bv = BV.size bv `seq` BV.nat bv `seq` ()
-}
-- instance NFData VertexType where rnf x = seq x ()

data VertexType = Root | Internal | Leaf | Network | Tree
    deriving (Read, Show, Eq) --NFData ?

type Vertex = (String, [Int], [Int], VertexType) -- name, child vertices, parent vertices, type (a bit redundant)
type Edge = (Int, Int, Maybe Double) -- terminal vertices (by numbers) and potential length
type PhyloGraphVect = (V.Vector Vertex, V.Vector Edge)
type GenPhyNetNode = (String, [String], [String]) --Make list for both so unresolved network (node name, [descendant name], [ancestor name])
type GenPhyNet = [GenPhyNetNode]

-- | null PhyloGraphVect
nullGraphVect :: PhyloGraphVect
nullGraphVect = (V.empty, V.empty) 

-- | 
-- Map a function over a traversable structure in parallel
-- Preferred over parMap which is limited to lists
-- Add chunking (with arguement) (via chunkList) "fmap blah blah `using` parListChunk chunkSize rseq/rpar"
-- but would have to do one for lists (with Chunk) and one for vectors  (splitAt recusively)
parmap :: Traversable t => Strategy b -> (a->b) -> t a -> t b
parmap strat f = withStrategy (parTraversable strat).fmap f

-- | seqParMap takes strategy,  if numThread == 1 retuns fmap otherwise parmap and 
seqParMap :: Traversable t => Strategy b -> (a -> b) -> t a -> t b
seqParMap strat f =
  if getNumThreads > 1 then parmap strat f
  else fmap f

-- | set myStrategy parallel strategy to rdeepseq
myStrategy :: (NFData b) => Strategy b
myStrategy = rdeepseq

-- | getNumThreads gets number of COncurrent  threads
{-# NOINLINE getNumThreads #-}
getNumThreads :: Int
getNumThreads = unsafePerformIO getNumCapabilities

-- | fgl2PGV takes an fgl (functional graph) and convertes to PhyloGraphVect
-- to use local (and old) Adams consensus functions
fgl2PGV :: P.Gr BV.BV (BV.BV, BV.BV) -> PhyloGraphVect
fgl2PGV inFunctionalGraph =
  (V.empty, V.empty)

-- | pgv2FGL take a PhyloGraphVect and converts to an fgl graph
pgv2FGL :: PhyloGraphVect -> P.Gr BV.BV (BV.BV, BV.BV)
pgv2FGL inPhyloGraphVect  =
  G.empty

-- | isTree takes fgl graph and checks is conected, no self loops, songle root, no indegree
-- > 1 nodes
isTree :: P.Gr BV.BV (BV.BV, BV.BV) -> Bool
isTree inFunctionalGraph = 
  True

-- | getRootNamesFromGenPhyNet extracts non-leaf-non-root 
-- names from vertices in order found
getRootNamesFromGenPhyNet :: GenPhyNet  -> [String]
getRootNamesFromGenPhyNet inNet =
    if null inNet then []
    else
        let (name, desc, anc) = head inNet
            vertType = getVertType (length desc) (length anc)
        in
        if (vertType == Root) then
            name : getRootNamesFromGenPhyNet (tail inNet)
        else 
            getRootNamesFromGenPhyNet (tail inNet)

-- | getNonLeafNonRootNamesFromGenPhyNet extracts non-leaf-non-root 
-- names from vertices in order found
getNonLeafNonRootNamesFromGenPhyNet :: GenPhyNet  -> [String]
getNonLeafNonRootNamesFromGenPhyNet inNet =
    if null inNet then []
    else
        let (name, desc, anc) = head inNet
            vertType = getVertType (length desc) (length anc)
        in
        if (vertType /= Leaf) && (vertType /= Root) then
            name : getNonLeafNonRootNamesFromGenPhyNet (tail inNet)
        else 
            getNonLeafNonRootNamesFromGenPhyNet (tail inNet)

-- | getLeafNamesFromGenPhyNet extracts leaf names from vertices in order found
getLeafNamesFromGenPhyNet :: GenPhyNet  -> [String]
getLeafNamesFromGenPhyNet inNet =
    if null inNet then []
    else
        let (name, desc, anc) = head inNet
        in
        if (getVertType (length desc) (length anc)) == Leaf then
            name : getLeafNamesFromGenPhyNet (tail inNet)
        else 
            getLeafNamesFromGenPhyNet (tail inNet)

-- | getVertType takes list of desc and anc to determine type of vertex
getVertType :: Int -> Int -> VertexType
getVertType nDesc nAnc =
    if (nDesc == 0) && (nAnc == 0) then error "Isolated node"
    else if nAnc == 0 then Root
    else if nDesc == 0 then Leaf
    else if nAnc == 1 then Tree
    else if nAnc > 2 then Network
    else error ("Screwey node: indegree " ++ show nDesc ++ " outdegree " ++ show nAnc)
 
-- | getVertNum takes a list of vertex names and teh complete list and
-- returns a list of the indices (integers) of the names
getVertNum :: [String] -> [String] -> [Int]
getVertNum nameList vertexNameList =
    if null vertexNameList then []
    else 
        let firstVertName = head vertexNameList
            vertNum = elemIndex firstVertName nameList
        in
        if vertNum == Nothing then error ("Error in vertex name index: " ++ show firstVertName ++ " in " ++ show nameList)
        else 
            (fromJust vertNum) : getVertNum nameList (tail vertexNameList)

-- | oldIndex2New takes a PhyloGraphVect and creates a list of reorder based ordered name list
oldIndex2New :: V.Vector Vertex  -> [String] -> [(Int, Vertex)]
oldIndex2New inVertexVect nameList =
    if V.null inVertexVect then []
    else
        let curVert = V.head inVertexVect
            (vertName, _, _, _) = curVert
            vertNum = elemIndex vertName nameList
        in
        (fromJust vertNum, curVert) : (oldIndex2New (V.tail inVertexVect) nameList)

-- | genForestToPhyloGraph converts GenForest to PhyloGraph (so can use legacy
-- ENewick etc parsers
-- takes flattened vector of GenPhyNetNodes and builds vertices (leaf, internal,
-- and root) and edges.  Vertices and edges are added to input null
-- PhyloGraphVect
-- EDGES SEEM TO BE INCORRECT IN PLACES
genForestToPhyloGraphVect :: V.Vector GenPhyNetNode -> PhyloGraphVect -> [String] -> PhyloGraphVect
genForestToPhyloGraphVect inGen inPhyVect nameList = 
    if V.null inGen then inPhyVect
    else 
        let (inVertexName, inVertexDescNameList, inVertexAncNameList) = V.head inGen
            (curVertVect, curEdgeVect) = inPhyVect
            descNumList = getVertNum nameList inVertexDescNameList 
            ancNumList = getVertNum nameList inVertexAncNameList
            vertType = getVertType (length descNumList) (length ancNumList)
            newEdgeVect = V.zip3 (V.fromList ancNumList) (V.replicate (length ancNumList) 
                (head $ getVertNum nameList [inVertexName])) 
                (V.replicate (length ancNumList) Nothing) --edge from anc to current, no weight info
        in
        genForestToPhyloGraphVect (V.tail inGen)
            (V.snoc curVertVect (inVertexName, descNumList, ancNumList, vertType), 
            curEdgeVect V.++ newEdgeVect) nameList 

-- | getNamesFromGenPhyNet extracts names from vertices in order found
-- leaves are first, then internal, root last
getNamesFromGenPhyNet :: GenPhyNet  -> [String]
getNamesFromGenPhyNet inNet =
    (sort $ getLeafNamesFromGenPhyNet inNet) ++ (getNonLeafNonRootNamesFromGenPhyNet inNet)
        ++ (getRootNamesFromGenPhyNet inNet)

-- | getShortestList takes list and length and list of lists and return
-- shortest list
getShortestList :: ([a], [b]) -> Int -> [([a],[b])] -> ([a],[b])
getShortestList bestList lengthBestList inListList =
    if null inListList then bestList
    else
        let curList = head inListList
            lengthCurList = length $ fst curList
        in
        if lengthCurList < lengthBestList then getShortestList curList lengthCurList (tail inListList)
        else  getShortestList bestList lengthBestList (tail inListList)

-- | getSplitList take an LUB, list of placed taxa, and vector of tree vertices
-- and returns a list of splits for each input tree (from tree vertices)
-- also filters out placed taxa
-- CLEANUP--many more ioperations than needed--should be passed as better
-- structure 
getSplitList :: [String] -> [String] -> ([[Int]], [[String]]) -> [[String]]
getSplitList curLUB placedTaxa (potentialVerts, vertLeafSet) = 
    if null curLUB then error "Null LUB in getSplitList"
    else 
       let  vertList =  [(x, y) | (x, y)  <- (zip vertLeafSet potentialVerts), intersect x curLUB == curLUB] 
            smallestVert = snd $ getShortestList (head vertList) (length $ fst $ head vertList) (tail vertList)
            vectVertLeafSet = V.fromList vertLeafSet --adds factor of "n", could pass another variable?
            rawLUBs = map (vectVertLeafSet V.!) smallestVert
            newLUBs = map (\\ placedTaxa) rawLUBs
        in
        newLUBs 

-- | replaceChar take set of charcters to be replaced by a char in a String
replaceChar :: [Char] -> Char -> Char -> Char
replaceChar inSet2Replace replaceChar2 inChar =
    if elem inChar inSet2Replace then replaceChar2
    else inChar

-- | getVertsFromIndexList takes a list of vertex vector indices and returns a list of
-- vertices
getVertsFromIndexList :: [Int] -> PhyloGraphVect -> [Vertex]
getVertsFromIndexList indexList inGraphVect  = 
    if null indexList then []
    else 
        let (vertVect, _) = inGraphVect
        in
        (vertVect V.! (head indexList)) :  (getVertsFromIndexList (tail indexList) inGraphVect)

-- | ggenPhyNet2PhyloGraphVect takes as input GenPhyNet and return 
-- PhyloGraphVect with root as last node
genPhyNet2PhyloGraphVect :: GenPhyNet -> PhyloGraphVect
genPhyNet2PhyloGraphVect inGenPhyNet =
    if null inGenPhyNet then error "Null GenPhyNet in genPhyNet2PhyloGraphVect"
    else 
        let nameList = getNamesFromGenPhyNet inGenPhyNet
            (vertVect, edgeVect) = genForestToPhyloGraphVect (V.fromList inGenPhyNet) nullGraphVect nameList
            newVertVect = vertVect V.// (oldIndex2New vertVect nameList)
        in
        (newVertVect, edgeVect)

-- | makeAdamsNodes takes root Adams node, rootLUB, vertex sets of input
-- trees and placed leaf set and constructs each Adams node in turn.
makeAdamsNodes :: GenPhyNet -> String -> [[String]] -> [String] -> [([[Int]], [[String]])] -> GenPhyNet
makeAdamsNodes inAdamsTree parentName inLUBList placedTaxa bothLeafLists = --inTreeVertexLists vertexLeafSetList =
   if null inLUBList then inAdamsTree
   else 
        let curLUB = head inLUBList
        in
        if length curLUB == 1 then --make nodes since done
            let newNode = (head curLUB, [], [parentName])
            in
            makeAdamsNodes (newNode : inAdamsTree) parentName (tail inLUBList) 
                (head curLUB : placedTaxa) bothLeafLists --inTreeVertexLists vertexLeafSetList
        else if length curLUB == 2 then
            let leftChild  = lub2TreeRep  [head curLUB]
                rightChild = lub2TreeRep  [last curLUB]
                newNode1 = (lub2TreeRep curLUB, [leftChild, rightChild], [parentName]) 
                newNode2 = (leftChild, [], [lub2TreeRep curLUB])
                newNode3 = (rightChild, [], [lub2TreeRep curLUB])
                newGenPhyNet = newNode2 : (newNode3 : (newNode1 : inAdamsTree))
                newPlacedTaxa = (lub2TreeRep curLUB) : (leftChild : (rightChild : placedTaxa))
            in
            makeAdamsNodes newGenPhyNet parentName (tail inLUBList)
                newPlacedTaxa bothLeafLists --inTreeVertexLists vertexLeafSetList
        else --core case with LUB creation and taxon placementg
            let splitListList = map (getSplitList curLUB placedTaxa) bothLeafLists --(zip inTreeVertexLists vertexLeafSetList) 
                newLUBpre = leastUpperBound splitListList
                newLUB =  map sort $ [x | x <- newLUBpre, length x > 0] --had "map sort $" was this "sort" necessary?  for List intersection?
                newNode = (lub2TreeRep curLUB, map lub2TreeRep newLUB, [parentName])
            in
            --trace ("New LUBs " ++ show newLUB ++ " newNode " ++ show newNode)
            (makeAdamsNodes (newNode : inAdamsTree) (lub2TreeRep curLUB) (tail inLUBList)
                placedTaxa bothLeafLists) ++ --inTreeVertexLists vertexLeafSetList) ++ 
                (makeAdamsNodes [] (lub2TreeRep curLUB) newLUB placedTaxa bothLeafLists) --inTreeVertexLists vertexLeafSetList)

-- | getLeafSetFromNodeName takes String name of node and returns sorted list of leaf
-- names--ASSUMES node names are not given in input and are assigned as trees
-- are parsed
getLeafSetFromNodeName :: Vertex -> [String]
getLeafSetFromNodeName inVertex =
    let (nodeName, _, _, _) = inVertex
    in
    if null nodeName then error "Null node name in getLeafSetFromNodeName"
    else 
        let rawList = map (replaceChar ['(', ')', ','] (' ')) nodeName
        in
        sort $ words rawList --this sort required

-- | lub2TreeRep takes list of names and makes into unresolved subtree 
-- in parens
lub2TreeRep :: [String] -> String
lub2TreeRep inStringList =
    if null inStringList then error "Null input in lub2TreeRep"
    else
        if length inStringList == 1 then head inStringList
        else
            let inside = init $ concat $ map  (++ ",") inStringList
            in
            ( '(' : inside ) ++ ")"

-- | getDecendantLeafList iputs a vertex and returns leaf set (as list of
-- leaf names as strings) descdended from
-- that vertex, if a leaf, returns that leaf
getDecendantLeafList :: [Vertex] -> PhyloGraphVect -> [String]
getDecendantLeafList inVertexList inGraphVect =
    if null inVertexList then []
    else
        let (curVertName, descList, _, vertType) = head inVertexList 
            descVertList = getVertsFromIndexList  descList inGraphVect
        in
        if vertType == Leaf then
            curVertName : (getDecendantLeafList (tail inVertexList) inGraphVect)
        else  
            (getDecendantLeafList [head descVertList] inGraphVect ) 
                ++ (getDecendantLeafList (tail descVertList) inGraphVect) 
                ++ (getDecendantLeafList (tail inVertexList) inGraphVect)

-- | getSplitLeafList takes a node and returns a list of list of descendent leaves
getSplitLeafList :: [Int] -> PhyloGraphVect -> [[String]]
getSplitLeafList descList inGraphVect =
    if null descList then []
    else 
        let curDesc = head descList
            (vertexVect, _)  = inGraphVect
            curLeaves = getDecendantLeafList [vertexVect V.! curDesc] inGraphVect
        in curLeaves : getSplitLeafList (tail descList) inGraphVect

-- | getSplitLeafListList takes list of descenndents for PhyloGraphVect and
-- returns a list of descendant list for each split of each tree
getSplitLeafListList :: [[Int]] -> [PhyloGraphVect] -> [[[String]]]
getSplitLeafListList descListList inGraphVectList =
    if null descListList then []
    else if null inGraphVectList then error "Diff numbers of descdent lists and graphs"
    else 
        let curIntList = head descListList
            curGraphVectList = head inGraphVectList
        in 
        getSplitLeafList curIntList curGraphVectList :
            (getSplitLeafListList (tail descListList) (tail inGraphVectList))

-- | lub2 takes two lists of lists of names and generates the pairswise set of
-- intersections
lub2 :: [[String]] -> [[String]] -> [[String]]
lub2 s1 s2 =
    if null s1 then []
    else if null s2 then []
    else
        let intersectFirst = (intersect (head s1) (head s2)) : (lub2 [head s1] (tail s2))
        in
        intersectFirst ++ (lub2 (tail s1) s2)

-- | leastUpperBound takes list of list vertex leaf descendants (as Strings) 
-- and returns LUB of Adams II (1972) consensus
leastUpperBound :: [[[String]]] -> [[String]]
leastUpperBound inVertexListList =
    if length inVertexListList < 2 then
        error "Too few name lists in leastUpperBound"
    else if length inVertexListList == 2 then 
        let x = head inVertexListList
            y = last inVertexListList
        in
        lub2 x y
    else
        let x = head inVertexListList
            y = head $ tail inVertexListList
            z = tail $ tail inVertexListList
            t = lub2 x y
        in
        leastUpperBound (t : z)  

-- | get second retriueves 2nd element of 4
getSecond :: (a, b, c, d) -> b
getSecond inTuple  = 
    let (_, b2, _, _) = inTuple
    in
    b2

-- | leafSetFromVertexVect takes vector of veritces and returns set of leaf
-- names
leafSetFromVertexVect :: Set.Set String -> V.Vector Vertex -> Set.Set String
leafSetFromVertexVect inSet inVerts =
    if V.null inVerts then inSet
    else 
        let (curName, _, _, curType) = V.head inVerts
        in
        if curType == Leaf then 
            leafSetFromVertexVect (Set.insert curName inSet) (V.tail inVerts)
        else 
            leafSetFromVertexVect inSet (V.tail inVerts)

-- | getLeafSet tgake a list pf PhyloGraphVect and returns a pair with 
-- True if the leaf sets are identical, and a list of the sets
getLeafSet :: PhyloGraphVect -> Set.Set String
getLeafSet inGraphVect =
    let (inVerts, _) = inGraphVect
    in
    leafSetFromVertexVect Set.empty inVerts  


-- | setEqual checks for set equality by difference between union and
-- intersection is empty
setEqual :: Ord a => Set.Set a -> Set.Set a-> Bool
setEqual firstSet secondSet =
    let combinedElem = Set.union firstSet secondSet
        sameElem = Set.intersection firstSet secondSet
    in
    Set.empty == (Set.difference combinedElem sameElem)


-- | getAndCheckLeafSets take graphs and checks that leaf sets are identical
getAndCheckLeafSets :: [PhyloGraphVect] -> (Bool, [Set.Set String])
getAndCheckLeafSets inGraphs =
    if null inGraphs then error "Empty graph list in getAndCheckLeafSets"
    else 
        let leafSetList = map getLeafSet inGraphs
            firstSet = head leafSetList
            setDiffList = map (setEqual firstSet) (tail leafSetList)
            allEmpty = all (True ==) setDiffList
        in 
        (allEmpty, leafSetList)

-- | findRoot take PhyloGraphVect and return root index and Vertex
findRoot :: Int -> PhyloGraphVect -> (Int, Vertex)
findRoot index inGraph  =
    let (vertexVect, _) = inGraph
    in
    if index < V.length vertexVect then
        let (_, _, _, vertexType) = vertexVect V.! index
        in
        if vertexType == Root then 
            (index, vertexVect V.! index)
        else 
            findRoot (index + 1) inGraph
     else error "Index exceeeds vertex number in findRoot"

-- | getAdamsIIPair inputs 2 PhyloGraphVects and returns AdamsII consensus 
getAdamsIIPair ::  PhyloGraphVect -> PhyloGraphVect -> PhyloGraphVect
getAdamsIIPair inGraphVectA inGraphVectB =
        let inGraphVectList = [inGraphVectA, inGraphVectB]
            (sameLeafSet, leafSets) = getAndCheckLeafSets inGraphVectList
            curVertexSets = map fst inGraphVectList
            rootPairList = map (findRoot 0) inGraphVectList
            --rootIndexList = map fst rootPairList
            rootVertexList = map snd rootPairList
            rootSplits = map getSecond rootVertexList
            rootSplitLeafListList = getSplitLeafListList rootSplits inGraphVectList
            rootLUBPre = leastUpperBound rootSplitLeafListList
            rootLUB = map sort $ [x | x <- rootLUBPre, length x > 0] --need map sort $
            
            --create nodes based on LUBs
            leavesPlaced = concat $  [x | x <- rootLUB, length x < 3]
            rootNode = ("root", map lub2TreeRep rootLUB, [])
            vertexLeafSetList = map (map getLeafSetFromNodeName) (map V.toList curVertexSets)
            potentialVertexSets = map (map getSecond) (map V.toList curVertexSets)
        in
        if (sameLeafSet == False) then error ("Leaf sets of input graphs do not match"  ++ show leafSets)
        else
          --return processed when have all nodes
            let allAdamsNodes = makeAdamsNodes [rootNode] "root" rootLUB leavesPlaced (zip potentialVertexSets vertexLeafSetList) --curVertexSets vertexLeafSetList
            in
            genPhyNet2PhyloGraphVect allAdamsNodes
            
-- | makeAdamsII takes a list of fgl graphs, convertes them to PhyloGraphVect 
-- makes the Adamns consensus and then converts back to fgl for return to EUN code
makeAdamsII :: [P.Gr BV.BV (BV.BV, BV.BV)] ->  P.Gr BV.BV (BV.BV, BV.BV)
makeAdamsII inFGList =
  if null inFGList then G.empty
  else 
    let allTreesList = seqParMap myStrategy isTree inFGList
        allTrees = foldl1' (&&) allTreesList
    in
    if not allTrees then error ("Input graphs are not all trees: " ++ (show allTreesList))
    else 
      let inPGVList = fmap fgl2PGV inFGList -- paralle problem with NFData seqParMap myStrategy fgl2PGV inFGList
          adamsPGV = foldl1' getAdamsIIPair inPGVList
      in
      pgv2FGL adamsPGV
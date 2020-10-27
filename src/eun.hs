 {- |
Module      :  eun.hs
Description :  Progam to calcualte Edge Union Network ala Miyagi and Wheeler 2019
               input graphviz dot files and newick
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

{-# LANGUAGE ScopedTypeVariables #-}

module Main where

import qualified Adams                             as A
import           Control.DeepSeq
import           Control.Parallel.Strategies
import qualified Data.Bits                         as B
import qualified Data.BitVector                    as BV
import           Data.Char
import qualified Data.Graph.Inductive.Graph        as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import qualified Data.Graph.Inductive.Query.BFS    as BFS
import           Data.GraphViz                     as GV
import           Data.GraphViz.Attributes.Complete (Attribute (Label),
                                                    Label (..))
import           Data.GraphViz.Commands.IO
import           Data.GraphViz.Printing
import           Data.List
import qualified Data.Map.Strict                   as Map
import           Data.Maybe
import           Data.Monoid
import qualified Data.Text.Lazy                    as T
import qualified Data.Vector                       as V
import qualified ParseCommands                     as PC
import qualified PhyloParsers                      as PhyP
import           System.Environment
import           System.IO
import qualified Data.Set                          as S
-- import Debug.Trace


-- NFData instance for parmap/rdeepseq
instance NFData BV.BV where
  rnf bv = BV.size bv `seq` BV.nat bv `seq` ()

-- |
-- Map a function over a traversable structure in parallel
-- Preferred over parMap which is limited to lists
parmap :: Traversable t => Strategy b -> (a->b) -> t a -> t b
parmap strat f = withStrategy (parTraversable strat).fmap f

-- | turnOnOutZeroBit turns on the bit 'nleaves" signifying that
-- the node is outdegree 1
-- this so outdegree one nodes and their child have differnet bit sets
turnOnOutZeroBit :: BV.BV -> Int -> BV.BV
turnOnOutZeroBit inBitVect nLeaves = BV.or [B.bit nLeaves, inBitVect]

-- | turnOffOutZeroBit turns off the bit 'nleaves" signifying that
-- the node is outdegree /= 1
-- this so outdegree one nodes and their child have differnet bit sets
turnOffOutZeroBit :: BV.BV -> Int -> BV.BV
turnOffOutZeroBit inBitVect nLeaves = BV.extract (nLeaves - 1) 0 inBitVect

-- | setOutDegreeOneBit takes the number of children, number of leaves and bitvector and
-- sets nleaves bit to 1 if number of children == 1 otherwise clears/removes nleaves bit
setOutDegreeOneBit :: Int -> Int -> BV.BV -> BV.BV
setOutDegreeOneBit outDegree nLeaves inBitVect =
  if outDegree == 1 then turnOnOutZeroBit inBitVect nLeaves
  else turnOffOutZeroBit inBitVect nLeaves

-- | getRoot takes a greaph and list of nodes and returns vertex with indegree 0
-- so assumes a connected graph--with a single root--not a forest
getRoots :: P.Gr String String -> [G.Node] -> [Int]
getRoots inGraph nodeList =
  if null nodeList then []
  else
    let firstNode = head nodeList
    in
    if (G.indeg inGraph firstNode == 0) && (G.outdeg inGraph firstNode > 0) then firstNode : getRoots inGraph (tail nodeList)
    else getRoots inGraph (tail nodeList)

-- | getUnconmnectedNOdes takes a graph and list of nodes and returns vertex with indegree 0
-- and outdegeee 0
getUnConnectedNodes :: P.Gr String String -> Int -> [G.Node] -> [G.LNode BV.BV]
getUnConnectedNodes inGraph nLeaves nodeList =
  if null nodeList then []
  else
    let firstNode = head nodeList
        newNode =  (firstNode, B.bit firstNode)
    in
    if G.deg inGraph firstNode == 0 then
      newNode : getUnConnectedNodes inGraph nLeaves (tail nodeList)
    else getUnConnectedNodes inGraph nLeaves (tail nodeList)


-- | makeNodeFromChildren gets bit vectors as union of children in a post order traversal from leaves
makeNodeFromChildren :: P.Gr String String -> Int -> V.Vector (G.LNode BV.BV) -> Int -> [G.LNode BV.BV]
makeNodeFromChildren inGraph nLeaves leafNodes myVertex =
  if myVertex < nLeaves then [leafNodes V.! myVertex]
    else
      let myChildren = G.suc inGraph myVertex
          myChildrenNodes = fmap (makeNodeFromChildren inGraph nLeaves leafNodes) myChildren
          myBV = setOutDegreeOneBit (length myChildren) nLeaves $ BV.or $ fmap (snd . head) myChildrenNodes
      in
      (myVertex, myBV) : concat myChildrenNodes

-- | getNodesFromARoot follows nodes connected to a root.
-- can be fmapped over roots to hit all--should be ok if multiple hits on nodes
-- since all labeled by BV.BVs  need to fuse them if multiple roots to make sure nodes are consistent
-- and only one per root--should be ok for multikple jhists of nodes since BVs are from childre
-- just wasted work.  Should nub after to maeksure only unique (by BV) nodes in list at end
getNodesFromARoot :: P.Gr String String -> Int -> [G.LNode BV.BV] -> Int -> [G.LNode BV.BV]
getNodesFromARoot inGraph nLeaves leafNodes rootVertex =
  if  G.isEmpty inGraph then error "Input graph is empty in getLabelledNodes"
  else
    let rootChildVerts = G.suc inGraph rootVertex

        -- recurse to children since assume only leaves can be labbeled with BV.BVs
        -- fmap becasue could be > 2 (as in at root)
        rootChildNewNodes = fmap (makeNodeFromChildren inGraph nLeaves (V.fromList leafNodes)) rootChildVerts

        rootBV = setOutDegreeOneBit (length rootChildVerts) nLeaves $ BV.or $ fmap (snd . head) rootChildNewNodes
    in
    (rootVertex, rootBV) : concat rootChildNewNodes

-- | getLabelledNodes labels nodes with bit vectors union of subtree leaves via post order traversal
-- adds nodes to reDoneNodes as they are preocessed
-- reorder NOdes is n^2 should be figured out how to keep them in order more efficeintly
getLabelledNodes :: P.Gr String String -> Int -> [G.LNode BV.BV] -> [G.LNode BV.BV]
getLabelledNodes inGraph nLeaves leafNodes  =
  -- trace ("getLabbeled graph with " ++ (show $ G.noNodes inGraph) ++ " nodes in " ++ (showGraph inGraph)) (
  if  G.isEmpty inGraph then error "Input graph is empty in getLabelledNodes"
  else
    let rootVertexList = getRoots inGraph (G.nodes inGraph)
        htuList = nub $ concatMap (getNodesFromARoot inGraph nLeaves leafNodes) rootVertexList
    in
     -- this for adding in missing data
    let unConnectedNodeList = getUnConnectedNodes inGraph nLeaves (G.nodes inGraph)
    in
    reorderLNodes (htuList ++ unConnectedNodeList)  0


-- | findLNode takes an index and looks for node with that as vertex and retuirns that node
findLNode :: Int -> [G.LNode BV.BV] -> G.LNode BV.BV
findLNode vertex lNodeList =
  if null lNodeList then error ("Node " ++ show vertex ++ " not found")
  else
      let (a,b) = head lNodeList
      in
      if a == vertex then (a,b)
      else findLNode vertex (tail lNodeList)

-- | reorderLNodes takes a list of nodes and reorders and order based on node vertex number
-- n^2 ugh
reorderLNodes :: [G.LNode BV.BV]  -> Int -> [G.LNode BV.BV]
reorderLNodes inNodeList index
  | null inNodeList = []
  | index == length inNodeList = []
  | otherwise =
    let newNode =  findLNode index inNodeList
    in
    newNode : reorderLNodes inNodeList (index + 1)
   -- )

-- | relabelEdgs creates (BV.BV, BV.BV) labnels for an edges
relabelEdge :: V.Vector (G.LNode BV.BV) -> G.LEdge String -> G.LEdge (BV.BV, BV.BV)
relabelEdge allNodesVect inLEdge =
  let (e,u,_) = inLEdge
      eNodeBV = snd (allNodesVect V.! e)
      uNodeBV = snd (allNodesVect V.! u)
  in
  (e,u,(eNodeBV,uNodeBV))

-- | getLeafNumber take Graph and gets nu,ber of leaves (outdegree = 0)
getLeafNumber :: P.Gr BV.BV (BV.BV, BV.BV) -> Int
getLeafNumber inGraph =
  let degOutList = G.outdeg inGraph <$> G.nodes inGraph
  in length $ filter (==0) degOutList

-- | findStrLabel checks Attributes (list f Attribute) from Graphvz to extract the String label of node
-- returns Maybe Text
findStrLabel :: Attributes -> Maybe T.Text
findStrLabel = getFirst . foldMap getStrLabel

-- | getStrLabel takes an Attribute and reurns Text if StrLabel found, mempty otherwise
getStrLabel :: Attribute -> First T.Text
getStrLabel (Label (StrLabel txt)) = First . Just $ txt
getStrLabel _                      = mempty

-- | getLeafString takes a pairs (node vertex number, graphViz Attributes)
-- and returns String name of leaf of Stringified nude number if unlabbeled
getLeafString :: (Int, Attributes) -> String
getLeafString (nodeIndex, nodeLabel) =
  let maybeTextLabel = findStrLabel nodeLabel
  in
  maybe (show nodeIndex) T.unpack maybeTextLabel

-- | getLeafList returns leaf complement of graph from DOT file
getLeafList ::  P.Gr Attributes Attributes -> [G.LNode String]
getLeafList inGraph =
  if G.isEmpty inGraph then []
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList
        (nodeVerts, _) = unzip leafList
        newLabels = fmap getLeafString leafList
        leafList' = zip nodeVerts newLabels
    in
    leafList'

-- | getLeafListNewick returns leaf complement of graph from newick file
-- difference from above is in the leaf label type
getLeafListNewick ::  P.Gr a b -> [G.LNode a]
getLeafListNewick inGraph =
  if G.isEmpty inGraph then []
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList
        (nodeVerts, _) = unzip leafList
        -- only different line
        newLabels = fmap snd leafList
        leafList' = zip nodeVerts newLabels
    in
    leafList'

-- | checkNodesSequential takes a list of nodes and returns booolean
-- True if nodes are input with seqeutial numerical indices
-- False if not--scres up reindexing later which assumes they are successive
checkNodesSequential :: G.Node -> [G.Node] -> Bool
checkNodesSequential prevNode inNodeList
  | null inNodeList = True
  | (head inNodeList - prevNode) /= 1 = False
  | otherwise = checkNodesSequential (head inNodeList) (tail inNodeList)

-- | reAnnotateGraphs takes parsed graph input and reformats for EUN
reAnnotateGraphs :: P.Gr String String -> P.Gr BV.BV (BV.BV, BV.BV)
reAnnotateGraphs inGraph =
  -- trace ("Reannotating " ++ (showGraph inGraph)) (
  if G.isEmpty inGraph then error "Input graph is empty in reAnnotateGraphs"
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        nLeaves = length $ filter (==0) degOutList
        leafVerts = [0..(nLeaves - 1)]
        leafIntegers = fmap B.bit leafVerts
        leafBitVects =  leafIntegers  -- fmap (BV.bitVec nLeaves) leafIntegers
        leafNodes = Prelude.zip leafVerts leafBitVects
        allNodes = getLabelledNodes inGraph nLeaves leafNodes
        allEdges = fmap (relabelEdge (V.fromList allNodes)) (G.labEdges inGraph)
    in
    -- assign HTU BV via postorder pass.
    G.mkGraph allNodes allEdges

-- | checkBVs looks at BV.BV of node and retuns FALSE if found True if not
checkBVs :: BV.BV -> [G.LNode BV.BV] -> Bool
checkBVs inBV nodeList =
  null nodeList || (
  let (_, bv) = head nodeList
  in
  inBV /= bv && checkBVs inBV (tail nodeList))

-- | checkBVs looks at BV.BV of node and retuns FALSE if found True if not
checkEdgeBVs :: (BV.BV, BV.BV) -> [G.LEdge (BV.BV, BV.BV)] -> Bool
checkEdgeBVs (inABV, inBBV) edgeList =
  null edgeList || (
  let (_, _, (aBV,bBV)) = head edgeList
  in
  not ((inABV == aBV) && (inBBV == bBV)) && checkEdgeBVs (inABV, inBBV) (tail edgeList))

-- | addAndReIndexUniqueNodes takes an inital list of nodes and adds new nodes reindexed
-- check identity by BV.BV
-- index bigger than size becasue starting after number of leaves
addAndReIndexUniqueNodes :: Int -> [G.LNode BV.BV] -> [G.LNode BV.BV] -> [G.LNode BV.BV]
addAndReIndexUniqueNodes newIndex nodesToExamine uniqueReIndexedNodes =
  if null nodesToExamine then uniqueReIndexedNodes
  else
    let (_, inBV) = head nodesToExamine
        isUnique = checkBVs inBV uniqueReIndexedNodes
    in
    if isUnique then
      let newNode = (newIndex, inBV)
      in
      addAndReIndexUniqueNodes (newIndex + 1) (tail nodesToExamine) (newNode : uniqueReIndexedNodes)
    else addAndReIndexUniqueNodes newIndex (tail nodesToExamine) uniqueReIndexedNodes


-- | getNodeIndex takes a BV.BV and returns a node with the same BV.BV
getNodeIndex :: BV.BV -> [G.LNode BV.BV] -> Int
getNodeIndex inBV nodeList =
  if null nodeList then error ("Node  with BV " ++ show inBV ++ " not found in getNodeIndex")
  else
    let (index, bv) = head nodeList
    in
    if bv == inBV then index
    else getNodeIndex inBV (tail nodeList)


-- | addAndReIndexEdges  takes list of indexed nodes and BV, a list of edges to examine and a list of edges to keep
-- checks for uniqueness of edges by BV.BVs on (e,u) and reindexes the edge nodes based on the node set with bit vectors
-- keep method either 'unique" or "all" to keep lists of unique or all edges
addAndReIndexEdges :: String -> [G.LNode BV.BV] -> [G.LEdge (BV.BV,BV.BV)] -> [G.LEdge (BV.BV,BV.BV)] -> [G.LEdge (BV.BV,BV.BV)]
addAndReIndexEdges keepMethod indexedNodes edgesToExamine uniqueReIndexedEdges =
  if null edgesToExamine then uniqueReIndexedEdges
  else
      let (_, _, (eBV, uBV)) = head edgesToExamine
          isUnique = checkEdgeBVs (eBV, uBV) uniqueReIndexedEdges
      in
      if (keepMethod == "all") || isUnique then
          -- Find nodes with BVs of edge
          let eNode = getNodeIndex eBV indexedNodes
              uNode = getNodeIndex uBV indexedNodes
              newEdge = (eNode, uNode, (eBV, uBV))
          in
          addAndReIndexEdges keepMethod indexedNodes (tail edgesToExamine) (newEdge : uniqueReIndexedEdges)
      else addAndReIndexEdges keepMethod indexedNodes (tail edgesToExamine) uniqueReIndexedEdges

-- | testEdge nodeList fullEdgeList) counter
testEdge :: [G.LNode BV.BV] -> [G.LEdge (BV.BV,BV.BV)] -> Int -> [G.LEdge (BV.BV,BV.BV)]
testEdge nodeList fullEdgeList counter =
  let firstPart = take counter fullEdgeList
      secondPart = drop counter fullEdgeList
      newEdgeList = firstPart ++ tail secondPart
      (e,u, l) = head secondPart
      (newGraph :: P.Gr BV.BV (BV.BV, BV.BV)) = G.mkGraph nodeList newEdgeList
      bfsNodes = BFS.bfs e newGraph
      foundU = find (== u) bfsNodes
  in
  [(e,u,l) | isNothing foundU]

-- | makeEUN take list of nodes and edges, deletes each edge (e,u) in turn makes graph,
-- checks for path between nodes e and u, if there is delete edge otherwise keep edeg in list for new graph
makeEUN ::  [G.LNode BV.BV] -> [G.LEdge (BV.BV,BV.BV)] -> P.Gr BV.BV (BV.BV, BV.BV)
makeEUN nodeList fullEdgeList =
  let counterList = [0..(length fullEdgeList - 1)]
      -- requiredEdges = concat $ fmap (testEdge nodeList fullEdgeList) counterList
      requiredEdges = concat $ parmap rdeepseq (testEdge nodeList fullEdgeList) counterList
      newGraph = G.mkGraph nodeList requiredEdges
  in
  newGraph

-- | getLeafLabelMatches tyakes the total list and looks for elements in the smaller local leaf set
-- retuns int index of the match or (-1) if not found so that leaf can be added in orginal order
getLeafLabelMatches ::[G.LNode String] -> G.LNode String -> (Int, Int)
getLeafLabelMatches localLeafList totNode =
  if null localLeafList then (-1, fst totNode)
  else
    let (index, leafString) = head localLeafList
    in
    if snd totNode == leafString then (index, fst totNode)
    else getLeafLabelMatches (tail localLeafList) totNode

-- | reIndexEdge takes an (Int, Int) map, labelled edge, and returns a new labelled edge with new e,u vertices
reIndexLEdge ::  Map.Map Int Int -> G.LEdge a -> G.LEdge String
reIndexLEdge vertexMap inEdge =
  if Map.null vertexMap then error "Null vertex map"
  else
    let (e,u,_) = inEdge
        newE = Map.lookup e vertexMap
        newU = Map.lookup u vertexMap
    in
    if isNothing newE then error ("Error looking up vertex " ++ show e ++ " in " ++ show (e,u))
    else if isNothing newU then error ("Error looking up vertex " ++ show u ++ " in " ++ show (e,u))
    else (fromJust newE, fromJust newU, "")


-- | reIndexAndAddLeaves takes rawGraphs and total input leaf sets and reindexes node, and edges, and adds
-- in leaves (with out edges) so that later processing can get bit vectors correct and match from
-- graph to graph.
-- new node set in teh total leaf set form all graphs plus teh local HTUs renumbered up based on added leaves
-- the map contains leaf mappings based on label of leaf, the HTUs extend that map with stright integers.
-- edges are re-indexed based on that map
reIndexAndAddLeavesEdges :: [G.LNode String] -> ([G.LNode String], P.Gr a a) -> P.Gr String String
reIndexAndAddLeavesEdges totallLeafSet (inputLeafList, inGraph) =
  if G.isEmpty inGraph then G.empty
  else
      -- reindex nodes and edges and add in new nodes (total leaf set + local HTUs)
      -- create a map between inputLeafSet and totalLeafSet which is the canonical enumeration
      -- then add in local HTU nodes and for map as well
      -- trace ("Original graph: " ++ (showGraph inGraph)) (
      let correspondanceList = fmap (getLeafLabelMatches inputLeafList) totallLeafSet
          matchList = filter ((/=(-1)).fst) correspondanceList
          --remove order dependancey
          -- htuList = [(length inputLeafList)..(length inputLeafList + htuNumber - 1)]
          htuList = fmap fst (G.labNodes inGraph) \\ fmap fst inputLeafList
          htuNumber =  length (G.labNodes inGraph) - length inputLeafList
          newHTUNumbers = [(length totallLeafSet)..(length totallLeafSet + htuNumber - 1)]
          htuMatchList = zip htuList newHTUNumbers
          vertexMap = Map.fromList (matchList ++ htuMatchList)
          reIndexedEdgeList = fmap (reIndexLEdge vertexMap) (G.labEdges inGraph)

          newNodeNumbers = [0..(length totallLeafSet + htuNumber - 1)]
          attributeList = replicate (length totallLeafSet + htuNumber) "" -- origAttribute
          newNodeList = zip newNodeNumbers attributeList
      in
      G.mkGraph newNodeList reIndexedEdgeList

-- | relabelNode takes nofde list and labels leaves with label and HTUs with String of HexCode of BV label
relabelNodes :: [G.LNode BV.BV] -> [G.LNode String] -> [G.LNode String]
relabelNodes inNodes leafLabelledNodes
  | null inNodes = []
  | not $ null leafLabelledNodes = head leafLabelledNodes : relabelNodes (tail inNodes) (tail leafLabelledNodes)
  | otherwise =
  let (vertex, _) = head inNodes
  in
  (vertex, "HTU" ++ show vertex) : relabelNodes (tail inNodes) []


-- | addGraphLabels take Graph and changes to add nodes labelled wiyth String, edges as well
addGraphLabels :: P.Gr BV.BV (BV.BV, BV.BV) -> [G.LNode String] -> P.Gr String String
addGraphLabels inGraph totallLeafSet
  | G.isEmpty inGraph = error "Empty graph in addGraphLabels"
  | null totallLeafSet = error "Empty leaf set in addGraphLabels"
  | otherwise =
  let newNodes = relabelNodes (G.labNodes inGraph) totallLeafSet
    -- newNodes = totallLeafSet ++ newHTUList
      (eList, uList) = unzip (G.edges inGraph)

      newEdges = zip3 eList uList (replicate (length eList) "")
  in
  -- trace ("Relabelled EUN : " ++ (showGraph $ G.mkGraph newNodes newEdges) ++ " from " ++ (show totallLeafSet))
  G.mkGraph newNodes newEdges

-- | getIntersectionEdges takes a node A and cretes directed edges to each other edge in [B]
-- with rulkesLEdge
--  if A intesect B = empty then no edge
--  else if A intesect B = B then create edge A->B
--  else if A intesect B = A then create edge B->A
--  else --can't happen
getIntersectionEdges ::[G.LNode BV.BV] -> G.LNode BV.BV -> [G.LEdge (BV.BV,BV.BV)]
getIntersectionEdges bNodeList aNode =
  if null bNodeList then []
  else
      let (aIndex, aBV) = aNode
          (bIndex, bBV) = head bNodeList
          intersection = BV.and [aBV, bBV]
      in
      if (intersection == 0) || (aBV == bBV) then getIntersectionEdges (tail bNodeList) aNode
      else (if intersection == aBV then (bIndex, aIndex, (bBV, aBV)) : getIntersectionEdges (tail bNodeList) aNode
      else if intersection == bBV then (aIndex, bIndex, (aBV, bBV)) : getIntersectionEdges (tail bNodeList) aNode
      else  getIntersectionEdges (tail bNodeList) aNode)

-- |  getThresholdNodes takes a threshold and keeps those unique objects present in the threshold percent or
-- higher.  Sorted by frequency (low to high)
getThresholdNodes :: Int -> Int -> [[G.LNode BV.BV]] -> [G.LNode BV.BV]
getThresholdNodes thresholdInt numLeaves objectListList
  | thresholdInt < 0 || thresholdInt > 100 = error "Threshold must be in range [0,100]"
  | null objectListList = error "Empty list of object lists in getThresholdObjects"
  | otherwise =
  let threshold = (fromIntegral thresholdInt / 100.0) :: Double
      numGraphs = fromIntegral $ length objectListList
      objectList = sort $ snd <$> concat objectListList
      objectGroupList = Data.List.group objectList
      indexList = [numLeaves..(numLeaves + length objectGroupList - 1)]
      uniqueList = zip indexList (fmap head objectGroupList)
      frequencyList = fmap (((/ numGraphs) . fromIntegral) . length) objectGroupList
      fullPairList = zip uniqueList frequencyList
  in
  -- trace ("There are " ++ (show $ length objectListList) ++ " to filter: " ++ (show uniqueList) ++ " " ++ (show frequencyList))
  fst <$> filter ((>= threshold). snd) fullPairList

-- |  getThresholdEdges takes a threshold and number of graphs and keeps those unique edges present in the threshold percent or
-- higher.  Sorted by frequency (low to high)
-- modified from getThresholdNodes due to type change in edges
-- used and number from numleaves so can use BV
getThresholdEdges :: (Show a, Ord a) => Int -> Int -> [a] -> [a]
getThresholdEdges thresholdInt numGraphsIn objectList
  | thresholdInt < 0 || thresholdInt > 100 = error "Threshold must be in range [0,100]"
  | null objectList = error "Empty list of object lists in getThresholdEdges"
  | otherwise =
  let threshold = (fromIntegral thresholdInt / 100.0) :: Double
      numGraphs = fromIntegral numGraphsIn
      objectGroupList = Data.List.group $ sort objectList
      uniqueList = fmap head objectGroupList
      frequencyList = fmap (((/ numGraphs) . fromIntegral) . length) objectGroupList
      fullPairList = zip uniqueList frequencyList
  in
  fst <$> filter ((>= threshold). snd) fullPairList

-- | getUnConnectedHTUs removes unconnected non-leaf nodes from graph
-- this could be done better by just taking teh vertecces in the used edges
-- that have a path to the leaf set?
getUnConnectedHTUs :: P.Gr a b ->  [G.LNode a] -> [G.Node]
getUnConnectedHTUs inGraph leafNodes
  | null leafNodes = error "Empty leaf node list in getConnectedHTUs"
  | G.isEmpty inGraph = error "Empty graph in getConnectedHTUs"
  | otherwise =
  let nLeaves = length leafNodes
      htuList = drop nLeaves $ G.nodes inGraph
      degOutList = fmap (G.deg inGraph) htuList
      newNodePair = zip degOutList htuList
      htuPairList = filter ((< 2).fst ) newNodePair
      (_, unConnectedLabHTUList) = unzip htuPairList
  in
  unConnectedLabHTUList


-- | getPostOrderVerts takes a vertex and traverses postorder to root places all visirted nodes in a set of found
-- vertices. Keeps placing new nodes in recursion list until a root is hit.  If a node is already in found set
-- it is not added to list of nodes to recurse 
-- returns set of visited nodes
getPostOrderVerts :: P.Gr BV.BV (BV.BV, BV.BV) -> S.Set G.Node -> [G.Node] -> S.Set G.Node
getPostOrderVerts inGraph foundVertSet inVertexList = 
  if null inVertexList then foundVertSet
  else 
    let firstVertex = head inVertexList
    in
    if S.member firstVertex foundVertSet then getPostOrderVerts inGraph foundVertSet (tail inVertexList)
    else 
      let newFoundSet = S.insert firstVertex foundVertSet 
          parentVerts = G.pre inGraph firstVertex
          -- parentIndexList = G.pre inGraph firstVertex
          -- parentLabelList = fmap fromJust $ fmap (G.lab inGraph) parentIndexList
          --parentVerts = zip parentIndexList parentLabelList
      in 
      getPostOrderVerts inGraph newFoundSet (inVertexList ++ parentVerts)

-- | verticesByPostorder takes a graph and a leaf set and an initially empty found vertex set
-- as the postorder pass takes place form each leaf, each visited vertex is placed in foundVertSet
-- when roots are hit, it recurses back untill all paths are traced to a root.
-- final final rgaph is created and retuyrned from foundVertSet and input list
-- could have edges unconnected to leaves if consistent edge leading to a subtree with inconsistent configuration
-- so are filtered out by making sure each vertex in an edge is in the vertex list
verticesByPostorder :: P.Gr BV.BV (BV.BV, BV.BV) -> [G.LNode BV.BV] ->  S.Set G.Node -> P.Gr BV.BV (BV.BV, BV.BV)
verticesByPostorder inGraph leafNodes foundVertSet =
  if G.isEmpty inGraph then error "Empty graph in verticesByPostorder"
  else if null leafNodes then 
    let vertexIndexList = S.toList foundVertSet
        vertexLabelList = fmap fromJust $ fmap (G.lab inGraph) vertexIndexList
        vertexList = zip vertexIndexList vertexLabelList
        edgeList = concat $ parmap rdeepseq (verifyEdge vertexIndexList) $ G.labEdges inGraph
    in G.mkGraph vertexList edgeList
  else 
    let firstLeaf = fst $ head leafNodes
        firstVertices = getPostOrderVerts inGraph foundVertSet [firstLeaf]
    in
    verticesByPostorder inGraph (tail leafNodes) (S.union foundVertSet firstVertices)

-- | verifyEdge takes a vertex index list and an edge and checks to see if 
-- the subtyending vertices are in the vertex list nad returns teh edge as asingleton list
-- if yes--else empty list (for mapping purposes)
verifyEdge :: [G.Node] -> G.LEdge (BV.BV, BV.BV) -> [G.LEdge (BV.BV, BV.BV)]
verifyEdge vertIndexList inEdge@(e,u,_) =
  if e `notElem` vertIndexList then []
  else if u `notElem` vertIndexList then []
  else [inEdge]

-- | removeUnconnectedHTUGraph iteratibvely removes unconnecvted HTU nodes iuntill graph stable
removeUnconnectedHTUGraph ::  P.Gr BV.BV (BV.BV, BV.BV) -> [G.LNode BV.BV] ->  P.Gr BV.BV (BV.BV, BV.BV)
removeUnconnectedHTUGraph inGraph leafNodes
  | G.isEmpty inGraph = error "Empty graph in removeUnconnectedHTUGraph"
  | null leafNodes = error "Empty leaf node list in removeUnconnectedHTUGraph"
  | otherwise =
  let newGraph = G.delNodes (getUnConnectedHTUs inGraph leafNodes) inGraph
  in
  if G.equal inGraph newGraph then inGraph
  else removeUnconnectedHTUGraph newGraph leafNodes


-- | sortInputArgs takes a list of arguments (Strings) nd retuns a pair of lists
-- of strings that are newick or graphviz dotFile filenames for later parsing
sortInputArgs :: [String] -> [String] -> ([T.Text],[T.Text],[String],[String],[String]) -> ([T.Text],[T.Text],[String],[String],[String])
sortInputArgs inContents inArgs (curFEN, curNewick, curDot, curNewFiles, curFENFILES) =
  if null inArgs then (curFEN, curNewick, curDot, curNewFiles, curFENFILES)
  else
    let firstFileName = head inArgs
        firstContents = filter (not . isSpace) $ head inContents
    in
    if head firstContents == '(' then
      sortInputArgs (tail inContents) (tail inArgs) (curFEN, T.pack firstContents : curNewick, curDot, firstFileName : curNewFiles, curFENFILES)
    else if head firstContents == '<' then
      sortInputArgs (tail inContents) (tail inArgs) (T.pack firstContents : curFEN, curNewick, curDot, curNewFiles, firstFileName : curFENFILES)
    else -- assumes DOT
      sortInputArgs (tail inContents) (tail inArgs) (curFEN, curNewick, firstFileName : curDot, curNewFiles, curFENFILES)

-- | nodeText2String takes a node with text label and returns a node with String label
nodeText2String :: G.LNode T.Text -> G.LNode String
nodeText2String (index, label) = (index, T.unpack label)

-- | fglTextA2TextString converts the graph types from Text A to Text String
fglTextB2Text :: (Show b) => P.Gr b Double -> P.Gr b T.Text
fglTextB2Text inGraph =
  if G.isEmpty inGraph then G.empty
  else
    let labNodes = G.labNodes inGraph
        labEdges = G.labEdges inGraph
        (eList, uList, labelList) = unzip3 labEdges
        --- newLabels = fmap toShortest labelList
        newLabels = fmap (T.pack . show) labelList
        newEdges = zip3 eList uList newLabels
    in
    G.mkGraph labNodes newEdges
  
-- | main driver
main :: IO ()
main =
  do
    -- Process arguments
    --  csv file first line taxon/leaf names, subsequent lines are distances--must be symmetrical
    args <- getArgs

    -- process args
    -- removed output format since ouotputting both "dot" and "FENewick" files
    let (method, threshold, _, outputFile, inputFileList) = PC.processCommands args

    -- let method = head args
    -- let threshold =  (read (args !! 1) :: Int)

    hPutStrLn stderr ("\nGraph combination method: " ++ method ++ " at threshold " ++ show threshold ++ "\n")

    -- Check for input Newick files and parse to fgl P.Gr graphs
    fileContentsList <- mapM readFile inputFileList -- (drop 2 args)
    let (forestEnhancedNewickFileTexts, newickFileTexts, dotArgs, newickArgs, forestEnhancedNewickArgs) = sortInputArgs fileContentsList inputFileList ([],[],[],[],[]) -- (drop 2 args) ([],[],[],[],[])
    hPutStrLn stderr ("Graphviz Files " ++ show dotArgs)
    hPutStrLn stderr ("Newick Files " ++ show newickArgs)
    hPutStrLn stderr ("Forest Enhanced Newick Files " ++ show forestEnhancedNewickArgs)

    --New Newick parser
    let newickGraphList = PhyP.forestEnhancedNewickStringList2FGLList (T.concat $ newickFileTexts ++ forestEnhancedNewickFileTexts)

    -- input dot files
    (dotGraphList :: [DotGraph G.Node]) <- mapM readDotFile dotArgs

    -- initial conversion to fgl Graph format
    let (inputGraphListDot :: [P.Gr Attributes Attributes]) = fmap dotToGraph  dotGraphList

    if (length dotGraphList + length newickGraphList) < 2 then error ("Need 2 or more input graphs and " ++ show (length dotGraphList + length newickGraphList) ++ " have been input")
    else hPutStrLn stderr ("\nThere are " ++ show (length dotGraphList + length newickGraphList) ++ " input graphs")
    -- Get leaf sets for each graph (dot and newick separately due to types) and then their union
    -- this to allow for missing data (leaves) in input graphs
    -- and leaves not in same order
    let inputLeafListsDot = parmap rdeepseq  getLeafList inputGraphListDot
    let inputLeafListsNewick = fmap nodeText2String <$> parmap rdeepseq  getLeafListNewick newickGraphList

    let totallLeafString = foldl' union [] (fmap (fmap snd) (inputLeafListsDot ++ inputLeafListsNewick))
    let totallLeafSet = zip [0..(length totallLeafString - 1)] totallLeafString
    hPutStrLn stderr ("There are " ++ show (length totallLeafSet) ++ " unique leaves in input graphs")

    let sanityListDot = parmap rdeepseq  (checkNodesSequential  (-1)) (fmap G.nodes inputGraphListDot)
    let sanityListNewick = parmap rdeepseq  (checkNodesSequential  (-1)) (fmap G.nodes newickGraphList)
    let sanityList = sanityListDot ++ sanityListNewick
    let allOK = foldl' (&&) True sanityList
    if allOK then hPutStrLn stderr "\nInput Graphs passed sanity checks"
    else error ("Sanity check error(s) on input graphs (False = Failed) Non-sequential node indices: " ++ show (zip [0..(length sanityList - 1)] sanityList))

    -- Add in "missing" leaves from individual graphs and renumber edges
    let fullLeafSetGraphsDot = parmap rdeepseq (reIndexAndAddLeavesEdges totallLeafSet) $ zip inputLeafListsDot inputGraphListDot
    let fullLeafSetGraphsNewick = parmap rdeepseq (reIndexAndAddLeavesEdges totallLeafSet) $ zip inputLeafListsNewick (fmap fglTextB2Text newickGraphList)
    let fullLeafSetGraphs = fullLeafSetGraphsDot ++ fullLeafSetGraphsNewick

    -- Reformat graphs with appropriate annotations, BV.BVs, etc
    let processedGraphs = parmap rdeepseq reAnnotateGraphs fullLeafSetGraphs -- inputGraphList

    -- Create lists of reindexed unique nodes and edges, identity by BV.BVs
    -- The drops to not reexamine leaves repeatedly
    -- Assumes leaves are first in list
    let numLeaves = getLeafNumber (head processedGraphs)
    let leafNodes = take numLeaves (G.labNodes $ head processedGraphs)
    let firstNodes = G.labNodes $ head processedGraphs
    let numFirstNodes = length firstNodes
    let unionNodes = sort $ leafNodes ++ addAndReIndexUniqueNodes numFirstNodes (concatMap (drop numLeaves) (G.labNodes <$> tail processedGraphs)) (drop numLeaves firstNodes)
    let unionEdges = addAndReIndexEdges "unique" unionNodes (concatMap G.labEdges (tail processedGraphs)) (G.labEdges $ head processedGraphs)

    -- Won't keep this--just for comparative purposes
    let (totalGraph :: P.Gr BV.BV (BV.BV, BV.BV)) = G.mkGraph unionNodes unionEdges
    hPutStrLn stderr ("\nTotal Graph with " ++ show (length $ G.labNodes totalGraph) ++ " nodes and " ++ show (length $ G.labEdges totalGraph) ++ " edges")

    -- EUN as in Miyagi anmd WHeeler (2019)
    let eunGraph = makeEUN unionNodes unionEdges
    let eunInfo =  "EUN deleted " ++ show (length unionEdges - length (G.labEdges eunGraph) ) ++ " of " ++ show (length unionEdges) ++ " total edges"
    -- add back labels for vertices and "GV.quickParams" for G.Gr String Double or whatever
    let labelledEUNGraph = addGraphLabels eunGraph totallLeafSet
    -- Create EUN Dot file
    let eunOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledEUNGraph -- eunGraph
    let eunOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph labelledEUNGraph] False

    -- Create strict consensus
    let intersectionBVs = foldl1' intersect (fmap (fmap snd . G.labNodes) processedGraphs)
    let numberList = [0..(length intersectionBVs - 1)]
    let intersectionNodes = zip numberList intersectionBVs
    let strictConInfo =  "There are " ++ show (length intersectionNodes) ++ " nodes present in all input graphs"
    let intersectionEdges = nub $ concatMap (getIntersectionEdges intersectionNodes) intersectionNodes
    let strictConsensusGraph = makeEUN intersectionNodes intersectionEdges
    let labelledConsensusGraph = addGraphLabels strictConsensusGraph totallLeafSet
    let strictConsensusOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledConsensusGraph
    let strictConsensusOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph labelledConsensusGraph] False

    -- Creat Adams II consensus
    let adamsII = A.makeAdamsII totallLeafSet fullLeafSetGraphs
    -- let labelledAdamsII = addGraphLabels adamsIIBV totallLeafSet
    let adamsIIInfo = "There are " ++ show (length $ G.nodes adamsII) ++ " nodes present in Adams II consensus"
    let adamsIIOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams adamsII
    let adamsIIOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph adamsII] False


    -- Create thresholdMajority rule Consensus and dot string
    -- vertex-based CUN-> Majority rule ->Strict
    let thresholdNodes = leafNodes ++ getThresholdNodes threshold numLeaves (fmap (drop numLeaves . G.labNodes) processedGraphs)
    let thresholdEdges = nub $ concat $ parmap rdeepseq (getIntersectionEdges thresholdNodes) thresholdNodes
    let thresholdConsensusGraph = makeEUN thresholdNodes thresholdEdges
    let thresholdConInfo =  "There are " ++ show (length thresholdNodes) ++ " nodes present in " ++ (show threshold ++ "%") ++ " of input graphs"
    let labelledTresholdConsensusGraph = addGraphLabels thresholdConsensusGraph totallLeafSet
    let thresholdConsensusOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledTresholdConsensusGraph
    let thresholdConsensusOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph labelledTresholdConsensusGraph] False

    -- Create threshold EUN and dot string (union nodes from regular EUN above)
    -- need to add filter to remove HTU degreee < 2
    -- need to add second pass EUN after adding HTU->Leaf edges
    let allEdges = addAndReIndexEdges "all" unionNodes (concatMap G.labEdges (tail processedGraphs)) (G.labEdges $ head processedGraphs)
    let thresholdEUNEdges = getThresholdEdges threshold (length processedGraphs) allEdges
    let thresholdEUNGraph' = makeEUN unionNodes thresholdEUNEdges

    -- Remove unnconnected HTU nodes via postorder pass from leaves
    let thresholdEUNGraph = verticesByPostorder thresholdEUNGraph' leafNodes S.empty
    -- let thresholdEUNGraph = removeUnconnectedHTUGraph thresholdEUNGraph' leafNodes
    let thresholdEUNInfo =  "\nThreshold EUN deleted " ++ show (length unionEdges - length (G.labEdges thresholdEUNGraph) ) ++ " of " ++ show (length unionEdges) ++ " total edges"
    -- add back labels for vertices and "GV.quickParams" for G.Gr String Double or whatever
    let thresholdLabelledEUNGraph = addGraphLabels thresholdEUNGraph totallLeafSet
    -- Create EUN Dot file
    let thresholdEUNOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams thresholdLabelledEUNGraph -- eunGraph
    let thresholdEUNOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph thresholdLabelledEUNGraph] False

    -- Output file
    let outDOT = outputFile ++ ".dot"
    let outFEN = outputFile ++ ".fen"
    if method == "eun" then
        if threshold == 0 then  do {hPutStrLn stderr eunInfo; writeFile outDOT eunOutDotString; writeFile outFEN eunOutFENString}
        else do {hPutStrLn stderr thresholdEUNInfo; writeFile outDOT thresholdEUNOutDotString; writeFile outFEN thresholdEUNOutFENString}
    else if method == "adams" then do {hPutStrLn stderr adamsIIInfo; writeFile outDOT adamsIIOutDotString; writeFile outFEN adamsIIOutFENString}
    else if method == "majority" then
        if threshold == 100 then do {hPutStrLn stderr strictConInfo; writeFile outDOT strictConsensusOutDotString; writeFile outFEN strictConsensusOutFENString}
        else do {hPutStrLn stderr thresholdConInfo; writeFile outDOT thresholdConsensusOutDotString; writeFile outFEN thresholdConsensusOutFENString}
    else if method == "strict" then do {hPutStrLn stderr strictConInfo; writeFile outDOT strictConsensusOutDotString; writeFile outFEN strictConsensusOutFENString}
    else error ("Graph combination method " ++ method ++ " is not implemented")

    hPutStrLn stderr "\nDone"



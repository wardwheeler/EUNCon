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

Todo:
  Add Adams consensus
  In/Output enewick?

-}

{-# LANGUAGE ScopedTypeVariables #-}

module Main where

import System.IO
import System.Environment
import Data.List
import Data.Maybe
import qualified Data.Vector as V
import Data.GraphViz as GV
import Data.GraphViz.Commands.IO
import Data.GraphViz.Printing
import Data.GraphViz.Attributes.Complete (Attribute (Label), Label (..))
import qualified Data.Graph.Inductive.Graph as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import qualified Data.Graph.Inductive.Query.BFS as BFS
import qualified Data.Text.Lazy as T
import qualified Data.Bits as B
import qualified Data.Map.Strict as Map
import Control.Parallel.Strategies
import Control.DeepSeq
import qualified Data.BitVector as BV
import Data.Monoid
import Data.Char
import qualified Adams as A
import qualified PhyloParsers as PhyP
-- import Data.Typeable

-- import Debug.Trace


-- NFData instance for parmap/rdeepseq
instance NFData BV.BV where
  rnf bv = BV.size bv `seq` BV.nat bv `seq` ()

-- | 
-- Map a function over a traversable structure in parallel
-- Preferred over parMap which is limited to lists
parmap :: Traversable t => Strategy b -> (a->b) -> t a -> t b
parmap strat f = withStrategy (parTraversable strat).fmap f


-- | functions for triples 
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

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


-- | writeFiles takes a stub and list of String
-- and writes to files usinf stub ++ number as naming convention
writeFiles :: String -> String -> Int -> [String] -> IO ()
writeFiles stub typeString number fileStuffList =
  if null fileStuffList then hPutStrLn stderr ("Wrote " ++ show number ++ " " ++ stub ++ ".X." ++ typeString ++ " files")
  else do
    writeFile (stub ++ "." ++ show number ++ "." ++ typeString) (head fileStuffList)
    writeFiles stub typeString (number + 1) (tail fileStuffList)

-- | getRoot takes a greaph and list of nodes and returns vertex with indegree 0
-- so assumes a connected graph--with a single root--not a forest
getRoots :: P.Gr String String -> [G.Node] -> [Int]
getRoots inGraph nodeList =
  if null nodeList then [] --error "Root vertex not found in getRoot"
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
        newNode =  (firstNode, B.bit firstNode) -- (firstNode, BV.bitVec nLeaves (2^firstNode))
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
          -- myBV = Main.or (fmap snd $ fmap head myChildrenNodes)
          -- newLNode = (myVertex, myBV)
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
  -- trace ((show inLEdge) ++ " " ++ (show (allNodesVect V.! e, allNodesVect V.! u)) ++ " " ++ "=>" ++ (show (e,u,(eNodeBV,uNodeBV))))
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


-- | getLeafList reutnds leaf complement of graph
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
getLeafListNewick ::  P.Gr String String -> [G.LNode String]
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
  -- )


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
    -- trace ("There are " ++ (show nLeaves) ++ " leaves ") -- ++ (show leaves))
    G.mkGraph allNodes allEdges
    -- )

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

-- | showGraph a semi-formatted show for Graphs
showGraph :: (Show a, Show b) => P.Gr a b -> String -- BV.BV (BV.BV, BV.BV) -> String
showGraph inGraph =
  if G.isEmpty inGraph then "Empty Graph"
  else
      let nodeString = show $ G.labNodes inGraph
          edgeString  = show $ G.labEdges inGraph
      in
      ("Nodes:" ++ nodeString ++ "\n" ++ "Edges: " ++ edgeString)

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

-- | isSelfEdge returns True if e=u in (e,u), else False
isSelfEdge :: G.LEdge a -> Bool
isSelfEdge (e,u,_) = e == u

-- | makeEdge takes a terminal vertex and its BV.BV and a parent vertex
-- and cretes a new labelled edge
makeEdge :: G.Node -> BV.BV ->  G.Node -> G.LEdge (BV.BV, BV.BV)
makeEdge uNode eBV eNode = (eNode, uNode, (eBV,eBV))


-- | makeOiutOneEdgeList takes a ;list of self edges and deg=0 nodes and makes new edges with nodes
-- the number should be equal
makeOutOneEdgeList :: P.Gr BV.BV (BV.BV,BV.BV) -> [G.LEdge (BV.BV,BV.BV)] -> [G.Node] -> [G.LEdge (BV.BV,BV.BV)]
makeOutOneEdgeList inGraph selfEdgeList unConnectedNodeList
  | length selfEdgeList > length unConnectedNodeList = error ("This seems like a problem--insufficent unconnected nodes: " ++ show unConnectedNodeList ++ " self edges: " ++ show selfEdgeList)
  | null selfEdgeList = []
  | otherwise =
    let (e,_,(eBV,_)) = head selfEdgeList
        u = head unConnectedNodeList
        newChildEdge = (u,e,(eBV,eBV))
        parentOfSelfEdgeNode = filter (/= e ) (G.pre inGraph e)
    in
    if null parentOfSelfEdgeNode then error ("No parents for self edge node " ++ show e)
    else
      let newParentEdges = fmap (makeEdge u eBV) parentOfSelfEdgeNode
      in
      (newChildEdge : newParentEdges) ++ makeOutOneEdgeList inGraph (tail selfEdgeList) (tail unConnectedNodeList)

-- | deleteEdges femoves edges in first leisty frpom edges in second
deleteEdges :: [G.LEdge (BV.BV,BV.BV)] -> [G.LEdge (BV.BV,BV.BV)]-> [G.LEdge (BV.BV,BV.BV)]
deleteEdges edgesToDelete edgesToCheckList
  | null edgesToDelete = edgesToCheckList
  | null edgesToCheckList = []
  | otherwise =
    let firstEdge = head edgesToCheckList
        toKeep = notElem firstEdge edgesToDelete
    in
    if not toKeep then deleteEdges edgesToDelete (tail edgesToCheckList)
    else firstEdge : deleteEdges edgesToDelete (tail edgesToCheckList)

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

-- areListsEqual etst two lists for equality in all elements
areListsEqual :: (Eq a) => [a] -> [a] -> Bool
areListsEqual inA inB
  | length inA /= length inB = False
  | null inA = True
  | otherwise =
    let firstA = head inA
        firstB = head inB
    in
    firstA == firstB && areListsEqual (tail inA) (tail inB)

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


-- | bvValue takes two nodes and compare based on BV for use in groupBy in getThresholdNodes
bvValue :: G.LNode BV.BV -> G.LNode BV.BV -> Bool
bvValue a b =
  snd a == snd b

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
sortInputArgs :: [String] -> [String] -> ([T.Text],[String],[String]) -> ([T.Text],[String],[String])
sortInputArgs inContents inArgs (curNewick, curDot, curNewFiles) =
  if null inArgs then (curNewick, curDot, curNewFiles)
  else
    let firstFileName = head inArgs
        firstContents = filter (not . isSpace) $ head inContents
    in
    if head firstContents == '(' then
      sortInputArgs (tail inContents) (tail inArgs) (T.pack firstContents : curNewick, curDot, firstFileName : curNewFiles)
    else
      sortInputArgs (tail inContents) (tail inArgs) (curNewick, firstFileName : curDot, curNewFiles)

-- | removeBranchLengths from Text group
removeBranchLengths :: T.Text -> T.Text
removeBranchLengths inName
  | T.null inName = inName
  | T.last inName == ')' = inName
  | not (T.any (==':') inName) = inName
  | otherwise = T.reverse $ T.tail $ T.dropWhile (/=':') $ T.reverse inName
    -- )

-- | removeNewickComments take string and removes all "[...]"
removeNewickComments :: T.Text -> T.Text
removeNewickComments inString
  | T.null inString = T.empty
  | not (T.any (==']') inString) = inString
  | otherwise =
  let firstPart = T.takeWhile (/='[') inString
      secondPart = T.tail $ T.dropWhile (/=']') inString
  in
  T.append firstPart (removeNewickComments secondPart)

-- | removeComma take string and removes leading and terminal commas if present
removeComma :: T.Text -> T.Text
removeComma inString
  | T.null inString = T.empty
  | not (T.any (==',') inString) = inString
  | T.head inString == ',' = removeComma (T.tail inString)
  | T.last inString == ',' = removeComma (T.init inString)
  | otherwise = inString

-- | getLeafNodesAndEdges thisNode leaves 
getLeafNodesAndEdges :: G.LNode String -> Int -> [T.Text] -> [G.LNode String] -> [G.LEdge String] -> ([G.LNode String],[G.LEdge String])
getLeafNodesAndEdges parentNode numNodes leaves inNodes inEdges
  | null leaves = (inNodes, inEdges)
  | T.null (head leaves) = getLeafNodesAndEdges parentNode numNodes (tail leaves) inNodes inEdges
  | otherwise =
  let leafName = T.unpack $ removeBranchLengths (head leaves)
      thisNode = (numNodes, leafName)
      thisEdge = (fst parentNode, numNodes, "Edge(" ++ show (fst parentNode) ++ "," ++ show numNodes ++ ")")
  in
  getLeafNodesAndEdges parentNode (numNodes + 1) (tail leaves) (thisNode : inNodes) (thisEdge : inEdges)

-- | countParensAndSplit counts the left and right parens and if they are equal and > 0 splits Text
countParensAndSplit :: T.Text -> Int -> Int -> Int -> (T.Text, T.Text)
countParensAndSplit groupText leftParen rightParen pos' =
  let pos = fromIntegral pos'
  in
  if T.null groupText then error "Group split not found"
  else if pos > T.length groupText then (groupText, T.empty)
  else
    let first = T.index groupText (pos - 1)  --T.head groupText
    in
    if first == '(' then
      -- leading terminal
      if (leftParen == rightParen) && (pos > 1) then (T.take (pos - 1) groupText, T.drop (pos - 1) groupText)
      else countParensAndSplit groupText (leftParen + 1) rightParen (pos' + 1)
    else if first == ')' then
      if leftParen == (rightParen + 1) then (T.take pos groupText, T.drop pos groupText)
      else countParensAndSplit groupText leftParen (rightParen + 1) (pos' + 1)
    else countParensAndSplit groupText leftParen rightParen (pos' + 1)

-- | splitNewickGroups groupText
splitNewickGroups :: T.Text -> [T.Text]
splitNewickGroups groupText =
  --trace ("Splitting " ++ (T.unpack groupText)) (
  if T.null groupText then []
  else
    let (first, rest) = countParensAndSplit groupText (0 :: Int) (0 :: Int) (1 :: Int)
    in
    --trace ("->" ++ (T.unpack $ removeComma first) ++ "|" ++ (T.unpack $ removeComma rest)) 
    removeComma first : splitNewickGroups (removeComma rest)
    --)

-- | splitLeafNonLeafText take alist of Text and splits into leaf (no '(') and
-- non-leaf (with '(') lists
splitLeafNonLeafText :: [T.Text] -> [T.Text] -> [T.Text] ->  ([T.Text],[T.Text])
splitLeafNonLeafText inTextList leafList nonLeafList =
  if null inTextList then (leafList, nonLeafList)
  else
    let firstGroup = head inTextList
        hasParen = T.any (=='(') firstGroup
    in
    if hasParen then splitLeafNonLeafText (tail inTextList) leafList (firstGroup : nonLeafList)
    else splitLeafNonLeafText (tail inTextList) (firstGroup : leafList) nonLeafList



-- | newickToGraph takes text of newick description (no terminal ';')
-- ther shoudl be no spaces in the text--filtered earlier.
-- and recurives creates a nodes and edges in fgl library style
-- progresively consumes the text, each apren block representing a node->subtree 
-- operates on lists of Text so can recurse over subgroups of whatever size
-- EDGES WRONG FOR test-3.tre
-- ??? make edges and nodes after split before recurse
newickToGraph :: [G.LNode String] -> [G.LEdge String] -> G.LNode String -> [T.Text] -> P.Gr String String
newickToGraph nodeList edgeList parentNode inNewickTextList =
  -- trace ("Parsing " ++ (show inNewickTextList))  (
  if null inNewickTextList then
      -- remove edge to root and make graph
      G.mkGraph nodeList (filter ((> (-1)).fst3) edgeList)
  else
      let inNewickText = head inNewickTextList
      in
      if T.null inNewickText then G.empty
      -- is a leaf name only
      else if not (T.any (=='(') inNewickText) then
        let leaves = T.split (==',') inNewickText
            (newNodes, newEdges) = getLeafNodesAndEdges parentNode (length nodeList) leaves [] []
        in
        newickToGraph (nodeList ++ newNodes) (edgeList ++ newEdges) parentNode (tail inNewickTextList)
      -- should have '(' ')' and some charcters
      else if T.length inNewickText < 3 then error ("Improper newick tree component: " ++ T.unpack inNewickText)
      else
        --remove  labels (e.g. branch length as :XXX) on total group defined by '(...)' and remove outer parens
        let groupText = removeBranchLengths $ T.tail $ T.init $ removeBranchLengths inNewickText
            isSubTree = T.any (=='(') groupText
            thisNode = (length nodeList, "HTU" ++ show (length nodeList))
            thisEdge = (fst parentNode, length nodeList, "Edge(" ++ show (fst parentNode) ++ "," ++ show (length nodeList) ++ ")")
        in
        -- a set of leaves
        if not isSubTree then
          let leaves = T.split (==',') groupText
              (newNodes, newEdges) = getLeafNodesAndEdges thisNode (1 + length nodeList) leaves [] []
          in
          newickToGraph ((thisNode : nodeList) ++ newNodes) ((thisEdge : edgeList) ++ newEdges) parentNode (tail inNewickTextList)

        -- is a sub tree
        else
          -- find subtrees and recurse
          let subTrees = splitNewickGroups groupText
              -- spit into leaf and non-leaf descenedents
              (leafNodes, nonLeafNodes) = splitLeafNonLeafText subTrees [] []
              -- make any leaf nodes and edges
              (newNodes, newEdges) = getLeafNodesAndEdges thisNode (1 + length nodeList) leafNodes [] []
          in
          -- recurse the remainder
          newickToGraph ((thisNode : nodeList) ++ newNodes) ((thisEdge : edgeList) ++ newEdges) thisNode (nonLeafNodes ++ tail inNewickTextList)
          


-- | main driver
main :: IO ()
main =
  do
    -- Process arguments
    --  csv file first line taxon/leaf names, subsequent lines are distances--must be symmetrical
    args <- getArgs
    if length args < 4 then error "Need at least 4 arguments: method (eun or connsesus), a minimum percent representation (Integer), and at least two input GraphViz dot or Newick files"
    else hPutStrLn stderr ("\nReading " ++ show (length args - 2) ++ " input files Output to stdout")
    Prelude.mapM_ (hPutStrLn stderr) $ fmap show (zip [0..(length args-3)] (drop 2 args))

    -- nmethod and min percent arg
    let method = head args
    let threshold =  (read (args !! 1) :: Int)

    hPutStrLn stderr ("\nGraph combination method: " ++ method ++ " at threshold " ++ show threshold ++ "\n")

    -- Check for input Newick files and parse to fgl P.Gr graphs
    fileContentsList <- mapM readFile (drop 2 args)
    let (newickFileTexts, dotArgs, newickArgs) = sortInputArgs fileContentsList (drop 2 args) ([],[],[])
    hPutStrLn stderr ("Graphviz Files " ++ show dotArgs)
    hPutStrLn stderr ("Newick Files " ++ show newickArgs)

    let newickTextList = filter (not.T.null) $ T.split (==';') $ removeNewickComments (T.concat newickFileTexts)
    -- hPutStrLn stderr $ show newickTextList
    let newickGraphList = fmap (newickToGraph [] [] (-1, "") . (:[])) newickTextList
    -- Prelude.mapM_ (hPutStrLn stderr) (fmap showGraph newickGraphList)

    -- input dot files
    (dotGraphList :: [DotGraph G.Node]) <- mapM readDotFile dotArgs

    -- initial conversion to fgl Graph format
    let (inputGraphListDot :: [P.Gr Attributes Attributes]) = fmap dotToGraph  dotGraphList

    -- Get leaf sets for each graph (dot and newick separately due to types) and then their union
    -- this to allow for missing data (leaves) in input graphs
    -- and leaves not in same order
    let inputLeafListsDot = parmap rdeepseq  getLeafList inputGraphListDot
    let inputLeafListsNewick = parmap rdeepseq  getLeafListNewick newickGraphList

    let totallLeafString = foldl' union [] (fmap (fmap snd) (inputLeafListsDot ++ inputLeafListsNewick))
    let totallLeafSet = zip [0..(length totallLeafString - 1)] totallLeafString
    hPutStrLn stderr ("There are " ++ show (length totallLeafSet) ++ " unique leaves in input graphs")

    let sanityListDot = parmap rdeepseq  (checkNodesSequential  (-1)) (fmap G.nodes inputGraphListDot)
    let sanityListNewick = parmap rdeepseq  (checkNodesSequential  (-1)) (fmap G.nodes newickGraphList)
    let sanityList = sanityListDot ++ sanityListNewick
    let allOK = foldl' (&&) True sanityList
    if allOK then hPutStrLn stderr "Input Graphs passed sanity checks"
    else error ("Sanity check error(s) on input graphs (False = Failed) Non-sequential node indices: " ++ show (zip [0..(length sanityList - 1)] sanityList))


    -- Add in "missing" leaves from individual graphs and renumber edges
    hPutStrLn stderr "Reindexing Dots"
    let fullLeafSetGraphsDot = parmap rdeepseq (reIndexAndAddLeavesEdges totallLeafSet) $ zip inputLeafListsDot inputGraphListDot
    -- Prelude.mapM_ (hPutStrLn stderr) (fmap showGraph fullLeafSetGraphsDot)
    hPutStrLn stderr "Reindexing Newicks"
    let fullLeafSetGraphsNewick = parmap rdeepseq (reIndexAndAddLeavesEdges totallLeafSet) $ zip inputLeafListsNewick newickGraphList
    -- Prelude.mapM_ (hPutStrLn stderr) (fmap showGraph fullLeafSetGraphsNewick)
    hPutStrLn stderr "Joining Graphs"
    let fullLeafSetGraphs = fullLeafSetGraphsDot ++ fullLeafSetGraphsNewick
    -- Prelude.mapM_ (hPutStrLn stderr) (fmap showGraph fullLeafSetGraphs)

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
    hPutStrLn stderr ("Total Graph with " ++ show (length $ G.labNodes totalGraph) ++ " nodes and " ++ show (length $ G.labEdges totalGraph) ++ " edges")

    -- EUN as in Miyagi anmd WHeeler (2019)
    let eunGraph = makeEUN unionNodes unionEdges
    let eunInfo =  "EUN deleted " ++ show (length unionEdges - length (G.labEdges eunGraph) ) ++ " of " ++ show (length unionEdges) ++ " total edges"
    -- add back labels for vertices and "GV.quickParams" for G.Gr String Double or whatever
    let labelledEUNGraph = addGraphLabels eunGraph totallLeafSet
    -- Create EUN Dot file 
    let eunOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledEUNGraph -- eunGraph



    -- Create strict consensus
    let intersectionBVs = foldl1' intersect (fmap (fmap snd . G.labNodes) processedGraphs)
    let numberList = [0..(length intersectionBVs - 1)]
    let intersectionNodes = zip numberList intersectionBVs
    let strictConInfo =  "There are " ++ show (length intersectionNodes) ++ " nodes present in all input graphs"
    let intersectionEdges = nub $ concatMap (getIntersectionEdges intersectionNodes) intersectionNodes
    let strictConsensusGraph = makeEUN intersectionNodes intersectionEdges
    let labelledConsensusGraph = addGraphLabels strictConsensusGraph totallLeafSet
    let strictConsensusOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledConsensusGraph

    -- Creat Adams II consensus
    let adamsIIBV = A.makeAdamsII processedGraphs
    let labelledAdamsII = addGraphLabels adamsIIBV totallLeafSet
    let adamsIIInfo = "There are " ++ show (length $ G.nodes labelledAdamsII) ++ " nodes present in Adams II consensus"
    let adamsIIOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledAdamsII


    -- Create thresholdMajority rule Consensus and dot string
    let thresholdNodes = leafNodes ++ getThresholdNodes threshold numLeaves (fmap (drop numLeaves . G.labNodes) processedGraphs)
    let thresholdEdges = nub $ concat $ parmap rdeepseq (getIntersectionEdges thresholdNodes) thresholdNodes
    let thresholdConsensusGraph = makeEUN thresholdNodes thresholdEdges
    let thresholdConInfo =  "There are " ++ show (length thresholdNodes) ++ " nodes present in " ++ (show threshold ++ "%") ++ " of input graphs"
    let labelledTresholdConsensusGraph = addGraphLabels thresholdConsensusGraph totallLeafSet
    let thresholdConsensusOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams labelledTresholdConsensusGraph

    -- Create threshold EUN and dot string (union nodes from regular EUN above)
    -- need to add filter to remove HTU degreee < 2 
    -- need to add second pass EUN after adding HTU->Leaf edges
    let allEdges = addAndReIndexEdges "all" unionNodes (concatMap G.labEdges (tail processedGraphs)) (G.labEdges $ head processedGraphs)
    let thresholdEUNEdges = getThresholdEdges threshold (length processedGraphs) allEdges
    let thresholdEUNGraph' = makeEUN unionNodes thresholdEUNEdges

    -- Remove unnconnected HTU nodes
    let thresholdEUNGraph = removeUnconnectedHTUGraph thresholdEUNGraph' leafNodes
    let thresholdEUNInfo =  "Threshold EUN deleted " ++ show (length unionEdges - length (G.labEdges thresholdEUNGraph) ) ++ " of " ++ show (length unionEdges) ++ " total edges"
    -- add back labels for vertices and "GV.quickParams" for G.Gr String Double or whatever
    let thresholdLabelledEUNGraph = addGraphLabels thresholdEUNGraph totallLeafSet
    -- Create EUN Dot file 
    let thresholdEUNOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams thresholdLabelledEUNGraph -- eunGraph

    -- Output Dot file
    if method == "eun" then
        if threshold == 0 then  do {hPutStrLn stderr eunInfo; putStrLn eunOutDotString}
        else do {hPutStrLn stderr thresholdEUNInfo; putStrLn thresholdEUNOutDotString}
    else if method == "adams" then do {hPutStrLn stderr adamsIIInfo; putStrLn adamsIIOutDotString}
    else if method == "strict" then
        if threshold == 100 then do {hPutStrLn stderr strictConInfo; putStrLn strictConsensusOutDotString}
        else do {hPutStrLn stderr thresholdConInfo; putStrLn thresholdConsensusOutDotString}
    else error ("Graph combination method " ++ method ++ " is not implemented")
    hPutStrLn stderr "Done"



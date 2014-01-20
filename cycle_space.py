'''
 Copyright (C) 2014 - Juan Pablo Carbajal

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

    (1) Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer. 

    (2) Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.  
    
    (3)The name of the author may not be used to
    endorse or promote products derived from this software without
    specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
'''

# Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

import networkx as nx
import numpy as np

def chords (G):
  """Return a new graph that contains the edges that are the chords of G.

      The chords are all the edges that are not in a spanning three of G.

  Parameters
  ----------
  G : graph
     A NetworkX graph.

  Returns
  -------
  C : A new graph with the chords of G.
  T : The spanning tree from which C was calculated.

  """
  if G.is_directed ():
    if G.is_multigraph ():
      T = nx.minimum_spanning_tree (nx.MultiGraph (G))
    else:
      T = nx.minimum_spanning_tree (nx.Graph (G))
  else:
    T = nx.minimum_spanning_tree (G)

  C     = G.copy ()
  edges = T.edges_iter ()

  for e in edges:
    try:
      C.remove_edge (*e)
    except:
      C.remove_edge (*e[::-1])

  #deg = C.degree_iter ();
  #for d in deg:
  #  if d[1] == 0:
  #    C.remove_node (d[0])

  # Recreate T to get the same type as G
  T = G.copy ()
  if G.is_multigraph ():
    edges = C.edges_iter (keys=True)
  else:
    edges = C.edges_iter ()

  for e in edges:
    T.remove_edge (*e)

  return T,C

def cycle_space (T,C):
  """Return a list of cycle graphs representing the fundamental cycles of the
     spanning tree T extended by the chords edges of the graph C

  Parameters
  ----------
  T: a tree graph
     A NetworkX graph.
  C: graph representing chords for the tree T.
    A NetworkX (multidi)graph.

  Returns
  -------
  Z : a list of cycles

  """
  Z     = list ();
  edges = C.edges_iter ();

  for idx,e in enumerate (edges):
    if T.has_edge(*e) or T.has_edge(*e[::-1]):
      Z.append (list (e))
    else:
      T.add_edge (*e)
      Z.append (nx.cycle_basis (nx.Graph(T))[0])
      T.remove_edge (*e)

  return Z

def cycle_space_matrix (G,T,C):
  """Return a the matrix describing the fundamental cycles in G.
     If G is not oriented and arbitrary orientation is taken.

  Parameters
  ----------
  T: a tree graph of G
  C: graph representing chords for the tree T.

  Returns
  -------
  M : Matrix of the fundametal cycles

  """

  nrow = len(G.edges())
  ncol = len(C.edges())
  M    = np.zeros ([nrow,ncol],dtype=np.int)

  if G.is_multigraph ():
    Cedges_iter = C.edges_iter (keys=True)
    Gedges = G.edges (keys=True)
    Tedges = T.edges (keys=True)
  else:
    Cedges_iter = C.edges_iter ()
    Gedges = G.edges ()
    Tedges = T.edges ()

  for col, e in enumerate(Cedges_iter):
    row = Gedges.index (e)
    M[row,col]  = 1

    try:
      einT  = T.edges ().index(e[:2])
      # The edge is in the tree with the same orientation, hence invert it
      row2        = Gedges.index (Tedges[einT])
      M[row2,col] = -1
    except ValueError:
      try:
        ieinT = T.edges ().index(e[:2][::-1])
        # The edge is in the tree with the opposite orientation, hence leave it
        row2        = Gedges.index (Tedges[ieinT])
        M[row2,col] = 1
      except ValueError:

        tmp = T.edges()
        tmp.append (e[:2])
        cyc = nx.cycle_basis (nx.Graph(tmp))[0]
        cyc_e = []

        for idx,node in enumerate (cyc):
          if  idx < len (cyc)-1:
            cyc_e.append ((node, cyc[idx+1]))
          else:
            cyc_e.append ((node, cyc[0]))

        if e[:2][::-1] in cyc_e:
          # The chord is in the cycle with the opposite orientation
          # invert all edges of the cycle
          cyc_e = [i[::-1] for i in cyc_e]

        elif e[:2] not in cyc_e:
          raise NameError ('Something went wrong! The edge {} is not in cycle {}'.format(e,cyc))

        cyc_e.remove (e[:2])
        for ce in cyc_e:
          try:
            einT  = T.edges ().index(ce)
            # The edge is in the tree with the same orientation
            row2        = Gedges.index (Tedges[einT])
            M[row2,col] = 1
          except ValueError:
            ieinT = T.edges ().index(ce[::-1])
            # The edge is in the tree with the opposite orientation
            row2        = Gedges.index (Tedges[ieinT])
            M[row2,col] = -1

  return M

def test():
  cables = [(1,2),(1,2),(1,2),(3,1),(3,2)]
  # Testing undirected graph
  Gu     = nx.Graph (cables)
  Cu, Tu = chords (Gu)
  Mu     = cycle_space_matrix (Gu,Tu,Cu)
  Zu     = cycle_space (Tu,Cu)

  # Testing undirected multigraph
  Gm     = nx.MultiGraph (cables)
  Cm, Tm = chords (Gm)
  Mm     = cycle_space_matrix (Gm,Tm,Cm)
  Zm     = cycle_space (Tm,Cm)

  # Testing directed graph
  Gd     = nx.DiGraph (cables)
  Cd, Td = chords (Gd)
  Md     = cycle_space_matrix (Gd,Td,Cd)
  Zd     = cycle_space (Td,Cd)

  # Testing directed multigraph
  Gmd     = nx.MultiDiGraph (cables)
  Cmd, Tmd = chords (Gmd)
  Mmd     = cycle_space_matrix (Gmd,Tmd,Cmd)
  Zmd     = cycle_space (Tmd,Cmd)

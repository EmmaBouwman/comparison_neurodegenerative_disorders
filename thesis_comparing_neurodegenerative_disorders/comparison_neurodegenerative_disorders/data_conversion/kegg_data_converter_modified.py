import networkx as nx # Python network module
from Bio.KEGG import REST 
import defusedxml.ElementTree as ET # For parsing the KGML file
import matplotlib.pyplot as plt # For making a simple plot of the result
from typing import List
import requests


class KeggDataConverter:
    """
    Class for converting a kegg file so that it can immediately be merged in Cytoscape with 
    STRING networks with Uniprot IDs. The KGML file contains protein complexes which will be 
    split into single protein units and their Uniprot identifiers are retrieved while also 
    retaining the original network structure. 

    NOTE: the conversion of Uniprot IDs causes a relatively longer runtime.
    """

    def __init__(self, kegg_pathway_id: str) -> None:
        """
        Constructor of KeggDataConverter class.

        Parameters:
          kegg_pathway_id (str): The identifier for the kegg pathway that has to be converted,
            e.g. "hsa05010" for the Alzheimer Pathway.
        """

        self.kegg_pathway_id = kegg_pathway_id
        self.nodes_with_multiple_names_dict = {}
        self.nodes_with_one_name_list = []

        self.kegg_graph = nx.Graph(name="kegg_graph")

    def fetch_kegg_data(self) -> str:
        base_url = "https://rest.kegg.jp/get/"
        url = f"{base_url}{self.kegg_pathway_id}/kgml"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            raise ValueError(f"Failed to fetch KEGG data: {response.status_code}")

    def convert_kgml_to_graph(self) -> None:
        """
        Parses through the kegg_pathway file and returns a nx.Graph structure of the pathway.

        Parameters:
          kegg_pathway (str): Identifier of the kegg pathway that will be converted.
        Returns:
          kegg_graph (nx.graph): A NetworkX undirected graph of the kegg pathway.
        """

        kgml_data = self.fetch_kegg_data()

        nodes_list = []
        edges_list = []

        # Get kgml file with kegg pathway from the KEGG website.
        # kgml_file = REST.kegg_get(self.kegg_pathway_id, "kgml")

        # Get kgml file from data files
        # kgml_file = "../data/hsa05010.xml" # The hsa05010_test.xml file can be used for testing purposes
        tree = ET.fromstring(kgml_data)
        
        # Loop over entries in kgml file and keep track of the nodes that do/don't 
        # contain multiple names
        for entry in tree.findall("entry"):
            entry_id = entry.attrib.pop("id")
            names_list = entry.get("name").split(" ")
            if len(names_list) > 1:
                self.nodes_with_multiple_names_dict[entry_id] = names_list
            else:
                self.nodes_with_one_name_list.append(entry_id)
            nodes_list.append((entry_id, entry.attrib))

        # Loop over relations in the kgml file, which indicate the edges of the graph. 
        for relation in tree.findall("relation"):
            first_node = relation.attrib.pop("entry1")
            second_node = relation.attrib.pop("entry2")

            # If an edge has multiple subtypes, they are all appended to the edge attributes.
            for i, subtype in enumerate(relation):
                relation.attrib[f"subtype{i+1}"] = subtype.attrib["name"]
                relation.attrib[f"value{i+1}"] = subtype.attrib["value"]
            edges_list.append((first_node, second_node, relation.attrib))
        
        self.kegg_graph.add_nodes_from(nodes_list)
        self.kegg_graph.add_edges_from(edges_list)
    

    def parse_nodes_and_modify_graph_structure(self) -> None:
        """
        Parses through the nodes with multiple names and splits these up while
        giving them new identifiers. Adds the Uniprot IDs to all nodes. Reconnects
        the edges between the nodes and adds edges that interconnect all nodes that
        were split from a single node.
        """

        for node in self.nodes_with_multiple_names_dict:
            new_id_list = []
            node_names = self.nodes_with_multiple_names_dict[node]

            for i, kegg_id in enumerate(node_names):
                # Create new id and attributes dict for the kegg identifier
                new_id = f"{node}_{i+1}"
                copy_of_attrib_dict = self.kegg_graph.nodes[node]
                copy_of_attrib_dict["name"] = kegg_id
                new_id_list.append(new_id)

                self.kegg_graph.add_nodes_from([(new_id, copy_of_attrib_dict)])
                self.add_uniprot_to_node_attributes(new_id, kegg_id)
                self.connect_node_edges_to_new_id(node, new_id)

            self.interconnect_new_nodes(new_id_list)
        
        for node in self.nodes_with_one_name_list:
            kegg_id = self.kegg_graph.nodes[node]["name"]
            self.add_uniprot_to_node_attributes(node, kegg_id)

        self.kegg_graph.remove_nodes_from(self.nodes_with_multiple_names_dict.keys())


    def interconnect_new_nodes(self, new_id_list: List[str]) -> None:
        """Interconnects the nodes with edges given by new_id_list. Regards these connections
        as binding/association connections since they were in a complex according to KEGG.
        
        Parameters: 
          new_id_list: list of node identifiers that need to all be connected with each other.
        Returns:
          None
        """

        interconnected_edge_dict = {"type": "PPrel", "subtype1": "binding/association", "value1": "---"}
        for id in new_id_list:
            interconnected_edges_list = [(id, x, interconnected_edge_dict) for x in new_id_list if x != id]
            self.kegg_graph.add_edges_from(interconnected_edges_list)
            new_id_list = new_id_list[1:]
    

    def connect_node_edges_to_new_id(self, node_id: str, new_id: str) -> None:
        """
        Connects every edge to "new_id" in kegg.graph that connects to "node_id". 

        Parameters:
          node (str): the id of the node that the edges connect to. 
          new_id (str): the id in kegg.digraph that all edges need to get connected to. 
        Returns: 
          None
        """

        copied_graph = self.kegg_graph.copy()
        for current_node, adjacent_nodes_dict in copied_graph.adjacency():
            for adjacent_node, edge_attrib in adjacent_nodes_dict.items():
                if current_node == node_id:
                    self.kegg_graph.add_edges_from([(new_id, adjacent_node, edge_attrib)])
                elif adjacent_node == node_id:
                    self.kegg_graph.add_edges_from([(current_node, new_id, edge_attrib)])


    def add_uniprot_to_node_attributes(self, node: str, kegg_id: str) -> None:
        """Adds the uniprot identifiers from the KEGG ids to the node attributes in self.kegg_graph.nodes"""
        #
        # print(REST.kegg_conv("uniprot", kegg_id).read())
        # if kegg_id[:3] == "hsa":
        #     kegg_request = REST.kegg_conv("uniprot", kegg_id).read()
        #     uniprot_id = kegg_request.split("\t")[1].replace('up:', '').split('\n')[0]
        #
        #     self.kegg_graph.nodes[node]["uniprot"] = uniprot_id

        try:
            response = REST.kegg_conv("uniprot", kegg_id).read()
            if response:
            # Parse the response to extract the Uniprot ID
                if kegg_id[:3] == "hsa":
                    uniprot_id = response.split("\t")[1].replace('up:', '').split('\n')[0]
                    self.kegg_graph.nodes[node]["uniprot"] = uniprot_id
            # else:
            #     print("Unexpected KEGG ID format:", kegg_id)
            #     # Handle unexpected KEGG ID format
            print(response)

        except Exception as e:
            print("Error fetching or processing Uniprot data:", e)
            # if response.status_code == 400:
            #     print(kegg_id)
            # Handle the error gracefully


    def save_graph_to_file(self, filename: str = "kegg_graph") -> None:
        """Saves the current self.kegg_graph to a gml file."""
        graph_for_file = nx.convert_node_labels_to_integers(self.kegg_graph, first_label=1)
        nx.write_gml(graph_for_file, f"../data/{filename}.gml")


    def draw_graph(self, name: str = "output") -> None:
        """Draws a graph of the current self.kegg_graph with matplotlib and saves it as png file"""

        nx.draw(self.kegg_graph, node_size=400, pos=nx.circular_layout(self.kegg_graph),
                node_color='green', edge_color='black', with_labels=True, font_size=8)

        print(self.kegg_graph.nodes)

        plt.savefig(f"{name}.png")

if __name__ == "__main__":
    alzheimer_pathway = KeggDataConverter("hsa05017")
    alzheimer_pathway.convert_kgml_to_graph()
    alzheimer_pathway.parse_nodes_and_modify_graph_structure()
    alzheimer_pathway.save_graph_to_file("sca_hsa05017")
    # alzheimer_pathway.draw_graph("hsa05010_pathway_figure")


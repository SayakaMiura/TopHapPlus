import sys

import os

import pydot

import copy

import argparse

import random

import numpy

import datetime

from Bio import Phylo

from Bio import AlignIO

from Bio.Seq import Seq

from Bio.Align import MultipleSeqAlignment





parser = argparse.ArgumentParser(description="Mutation pathmaker.")

parser.add_argument("aln", help="One-hot encoded sequence alignment file.", type=str)

parser.add_argument("labels", help="One-hot encoded feature labels.", type=str)

parser.add_argument("--threshold", help="Mutation co-occurrence threshold.", type=float, default=0.75)

parser.add_argument("--flip_pass_thresh", help="Sets proportion of affirmative flip votes needed to pass (0.5 - 1.0).", type=float, default=0.66)

parser.add_argument("--flip_vote_thresh", help="Sets paradox score reduction needed for a pair to vote to flip (0.0 - 1.0).", type=float, default=0.9)

parser.add_argument("-o", "--output", help="Output directory to put results in.", type=str, default=".")

parser.add_argument("--all_times", help="List of all sampling times for complementing.", type=str, default=None)

parser.add_argument("--bootstrap_reps", help="Number of replicates to perform bootstrapping with.", type=int, default=0)

parser.add_argument("--bootstrap_ratio", help="Number of sequences to sample for each bootstrap replicate, specified as proportion of total number of input sequences.", type=float, default=1)

parser.add_argument("--stats_file", help="File to dump performance metrics to.", type=str, default=None)

parser.add_argument("--disable_graph_flipping", help="Turn off graph child-count-based mutation flipping.", action='store_true', default=False)

parser.add_argument("--initial_tree", help="Initial mutation tree to start with.", type=str, default=None)

parser.add_argument("--lock_init_tree", help="Prevent nodes that were in initial tree from being regrafted.", action='store_true', default=False)

parser.add_argument("--spr_mismatch_weight", help="Regrafting pass mismatch weight:\n"

												  "\t>1.0: Favor reducing mismatches over increasing matches."

												  "\t<1.0: Favor increasing matches over reducing mismatches", type=float, default=1.0)

vlevel = 2



args = parser.parse_args()





if args.output != "." and not os.path.exists(args.output):

	os.mkdir(args.output)



presence_char = "T"

absence_char = "A"

missing_data_char = "?"

verbose_stats = True

co_occurence_thresh = args.threshold

output_basename = os.path.splitext(os.path.basename(args.aln))[0]

cache_bootstrap_aln = True





def parse_initial_tree(tree_file):

	graph = pydot.graph_from_dot_file(tree_file)[0]

	comped_nodes = set()

	# Get all main edges from dot graph

	edges = [(x.get_source().replace("\"", ""), x.get_destination().replace("\"", "")) for x in graph.get_edges() if x.obj_dict['attributes'].get("constraint", "True") == "True"]

	for edge in edges:

		if "_comp" in edge[0]:

			comped_nodes.update([edge[0].replace("_comp","")])

		if "_comp" in edge[1]:

			comped_nodes.update([edge[1].replace("_comp", "")])

	edges = [(edge[0].replace("_comp",""), edge[1].replace("_comp","")) for edge in edges]

	return recur_children(edges, "root"), comped_nodes





def recur_children(edges, target):

	root = (target, [])

	for edge in edges:

		if edge[0] == target:

			root[1].append(recur_children(edges, edge[1]))

	return root





def build_comatrix_old(input_aln, labels):

	comatrix = {}

	comatrix_counts = {}

	for i in range(0, len(labels)):

		comatrix[labels[i]] = {}

		comatrix_counts[labels[i]] = {}

		for j in range(0, len(labels)):

			# print("{}//{}".format(labels[i],labels[j]))

			comatrix[labels[i]][labels[j]] = float(sum([1 if seq[i] + seq[j] == presence_char + presence_char else 0 for seq in input_aln])) / float(

				sum([1 if seq[j] != missing_data_char and seq[i] == presence_char else 0 for seq in input_aln]))

			comatrix_counts[labels[i]][labels[j]] = sum([1 if seq[i] + seq[j] == presence_char + presence_char else 0 for seq in input_aln])

	for key1 in comatrix.keys():

		comatrix[key1]["root"] = 1.0

		# for key2 in comatrix.keys():

		# 	comatrix[key1].update({key2: comatrix[key1][key2]})

		# 	comatrix[key2].update({key1: comatrix[key1][key2]})

		# 	comatrix_counts[key1].update({key2: comatrix_counts[key1][key2]})

		# 	comatrix_counts[key2].update({key1: comatrix_counts[key1][key2]})

	return comatrix, comatrix_counts





def build_comatrix(input_aln, labels):

	comatrix_dict = {}

	comatrix_counts_dict = {}

	presence_comatrix_list = []

	missing_comatrix_list = []

	for i in range(0, len(input_aln)):

		seq = str(input_aln[i].seq)

		bin_seq_present = [1 if char == presence_char else 0 for char in seq]

		bin_seq_missing = [0 if char == missing_data_char else 1 for char in seq]

		#seq_comatrix_present = numpy.outer(bin_seq_present, bin_seq_present)

		presence_comatrix_list.append(numpy.outer(bin_seq_present, bin_seq_present))

		missing_comatrix_list.append(numpy.outer(bin_seq_present, bin_seq_missing))

	presence_comatrix_counts = numpy.sum(presence_comatrix_list, 0)

	missing_comatrix_counts = numpy.sum(missing_comatrix_list, 0)

	for i in range(0,len(labels)):

		comatrix_dict[labels[i]] = {"root": 1.0}

		comatrix_counts_dict[labels[i]] = {}

		for j in range(0, len(labels)):

			if presence_comatrix_counts[i][j] < 1:

				comatrix_dict[labels[i]][labels[j]] = 0.0

			elif missing_comatrix_counts[i][j] < 1:

				raise Exception("Division by zero at indices: {}(i), {}(j)".format(i, j))

			else:

				comatrix_dict[labels[i]][labels[j]] = presence_comatrix_counts[i][j]/missing_comatrix_counts[i][j]

			comatrix_counts_dict[labels[i]][labels[j]] = presence_comatrix_counts[i][j]

	return comatrix_dict, comatrix_counts_dict





def find_lowest_parent(target, start_node, comatrix):

	for child in start_node[1]:

		#if comatrix[target][target] * co_occurence_thresh <= comatrix[target].get(child[0], comatrix[child[0]][target]):

		if co_occurence_thresh <= comatrix[target][child[0]]:

			return find_lowest_parent(target, child, comatrix)

	return start_node





def build_tree(labels, comatrix):

	# tree tuples have format (mut_pos,[child1, child2, ...])

	root = ("root", [])

	for i in range(0, len(labels)):

		parent = find_lowest_parent(labels[i], root, comatrix)

		if parent is None:

			root[1].append((labels[i], []))

		else:

			parent[1].append((labels[i], []))

	return root





def update_tree(labels, comatrix, root):

	existing_nodes = recur_nodes(root)

	for i in range(0, len(labels)):

		if labels[i] in existing_nodes:

			continue

		parent = find_lowest_parent(labels[i], root, comatrix)

		if parent is None:

			root[1].append((labels[i], []))

		else:

			parent[1].append((labels[i], []))

	return root





def recur_nodes(root):

	nodes = [root[0]]

	for subroot in root[1]:

		nodes.extend(recur_nodes(subroot))

	return nodes





def get_parent(root, target):

	parent = root

	children = root[1]

	found_parent = False

	if target in [x[0] for x in children]:

		found_parent = True

	while not found_parent:

		new_children = []

		for child in children:

			if target in [x[0] for x in child[1]]:

				parent = child

				found_parent = True

				break

			new_children.extend(child[1])

		children = new_children

	return parent





def get_descendants(node, descendants=None, level=0):

	if descendants is None:

		descendants = []

	for child in node[1]:

		descendants.append(child[0])

		get_descendants(child, descendants, level + 1)

	if level == 0:

		return descendants





def get_depth(root, target):

	depth = 1

	children = root[1]

	while True:

		if target in [child[0] for child in children]:

			break

		depth += 1

		new_children = []

		for child in children:

			new_children.extend(child[1])

		children = new_children

	return depth





def get_node_dict(root, node_dict=None, level=0):

	if node_dict is None:

		node_dict = {}

	node_dict.update({root[0]: root})

	for child in root[1]:

		get_node_dict(child, node_dict, level + 1)

	if level == 0:

		return node_dict





def parse_node_label(label_string):

	return {val.split(":")[0]: int(val.split(":")[1].strip()) for val in label_string.split("\nDropped:")[0].split('\n')[1:]}





def reset_highlighting(input_graph):

	for object in input_graph.get_nodes():

		if "color" in object.obj_dict["attributes"].keys():

			del object.obj_dict["attributes"]["penwidth"]

			del object.obj_dict["attributes"]["color"]

	for object in input_graph.get_edges():

		if "color" in object.obj_dict["attributes"].keys():

			if object.obj_dict['attributes']['color'] == "green":

				continue

			input_graph.del_edge(object.get_source(), object.get_destination())





def binary_complement_positions(initial_aln, pos_list):

	print("Binary complementing: {}".format(pos_list))

	for i in range(0, len(initial_aln)):

		seq_array = [initial_aln[i].seq[j] for j in range(0, len(initial_aln[i].seq))]

		for pos in pos_list:

			if seq_array[pos] == presence_char:

				seq_array[pos] = absence_char

			elif seq_array[pos] == absence_char:

				seq_array[pos] = presence_char

		initial_aln[i].seq = Seq("".join(seq_array) )#, initial_aln[i].seq.alphabet)





def regraft(graph, root, labels, target_node, new_parent):

	node_dict = get_node_dict(root)

	temp_clade = copy.deepcopy(node_dict[target_node])

	old_parent = get_parent(root, target_node)

	old_parent[1].remove(node_dict[target_node])

	node_dict[new_parent][1].append(temp_clade)

	if old_parent[0] != "root":

		old_src = labels.index(old_parent[0])

	else:

		old_src = "root"

	if vlevel >= 2:

		print("Removing edge ({}({}) -> {}({}))".format(old_parent[0], old_src, target_node, str(labels.index(target_node))))

	graph.del_edge(old_src, labels.index(target_node))

	if vlevel >= 2:

		print("Adding edge ({}({}) -> {}({}))".format(new_parent, str(labels.index(new_parent)), target_node, str(labels.index(target_node))))

	graph.add_edge(pydot.Edge(str(labels.index(new_parent)), labels.index(target_node), color="green", penwidth=3.0))





def graph_to_ancestors_dict(graph):

	edges = [(str(edge.get_source()), str(edge.get_destination())) for edge in graph.get_edges() if

				"color" not in edge.obj_dict['attributes'].keys() or edge.obj_dict['attributes']['color'] == "green"]

	node_labels = {node.obj_dict['name']: node.obj_dict['attributes']['label'].split('\n')[0].strip() for node in graph.get_nodes()}

	paths = {}

	trace_all_paths(edges, ["root"], paths)

	labeled_paths = {}

	for key in paths.keys():

		path = [node_labels[node] for node in paths[key]]

		labeled_paths[node_labels[key]] = path

	return labeled_paths





def trace_all_paths(edges, path, path_list=None):

	if path_list is None:

		path_list = {}

	for edge in edges:

		if edge[0] == path[-1]:

			trace_all_paths(edges, path + [edge[1]], path_list)

	path_list[copy.deepcopy(path[-1])] = copy.deepcopy(path)





class JobBase:

	def __init__(self, params):

		start_time = datetime.datetime.now()

		self.alignment = AlignIO.read(open(params.aln), "fasta")

		self.mut_labels = []

		self.ro_labels = []

		self.mut_days = {}

		self.params = params

		self.is_bootstrapping = False

		with open(params.labels, 'r') as file:

			for line in file:

				data = line.strip().split('\t')

				if data[0] == "label":

					continue

				self.mut_labels.append(data[0])

				self.mut_days[data[0]] = [int(x) for x in data[1].split(",")]

		if self.alignment.get_alignment_length() != len(self.mut_labels):

			print("Length of alignment not equal to number of mutation labels provided, stopping...")

			exit()

		self.species_count = len(self.alignment)

		self.aln_len = len(self.mut_labels)

		self.comped_nodes = set()

		if params.initial_tree is not None:
			
			(self.initial_tree, self.comped_nodes) = parse_initial_tree(params.initial_tree)

			if len(self.comped_nodes) > 0:

				binary_complement_positions(self.alignment, [self.mut_labels.index(node) for node in self.comped_nodes])

			self.initial_nodes = recur_nodes(self.initial_tree)

			missing_nodes = [node for node in self.initial_nodes if node not in self.mut_labels]

			missing_nodes.remove("root")
            #print (missing_nodes,self.initial_nodes,self.mut_labels)
			if len(missing_nodes) > 0:

				raise Exception("Missing data for the following nodes in initial tree:\n" + "\n".join(missing_nodes))
		#	print (missing_nodes,self.initial_nodes,self.mut_labels)
	#		open('a','r').readlines()
			# for node in recur_nodes(self.initial_tree):

			# 	print(node)

			# 	if node in self.comped_nodes:

			# 		binary_complement_positions()

		if vlevel > 1:

			print("Initial parsing time:{}".format(datetime.datetime.now() - start_time))

		self.ro_aln = self.alignment

		(self.main_results, self.modified_inputs) = self.main_analysis()

		rcount = 0

		for result in self.main_results:

			rcount += 1

			self.add_timing_labels(result, self.mut_days)

			result.write(os.path.join(self.params.output, output_basename + "_{}.txt".format(rcount)))

			result.write_png(os.path.join(self.params.output, output_basename + "_{}.png".format(rcount)))

			if rcount == len(self.main_results):

				edge_count = 0

				co_occurence_total = 0

				for edge in result.get_edges():

					if "label" in edge.obj_dict['attributes'].keys():

						edge_count += 1

						co_occurence_total += float(edge.obj_dict['attributes']['label'].strip().replace("%", ""))

				print("Average child/parent cooccurrence: {}% ({} edges)".format(round(co_occurence_total/edge_count, 1), edge_count))

				if self.params.stats_file is not None:

					with open(self.params.stats_file, 'a+') as file:

						file.write("{}\t{}\t{}\t".format(params.aln, round(co_occurence_total/edge_count, 1), edge_count))

		self.assign_seqs_to_nodes(self.modified_inputs[0], self.modified_inputs[1], self.main_results[-1])

		self.bs_results_list = []

		if self.params.bootstrap_reps > 0:

			self.is_bootstrapping = True

			self.species_count = int(self.params.bootstrap_ratio * len(self.alignment))

			self.bootstrap(self.params.bootstrap_reps, int(self.species_count))

			edge_counts = self.analyze_bootstrap_edges()

			self.main_results[-1].write(os.path.join(self.params.output, output_basename + "_{}.png".format(rcount)))

			self.main_results[-1].write_png(os.path.join(self.params.output, output_basename + "_{}.png".format(rcount)))

			with open(os.path.join(self.params.output, output_basename + "_{}_bootstrap_edge_counts.txt".format(len(self.main_results))), 'w') as file:

				edges = edge_counts.keys()

				nodes = list(set([edge[0] for edge in edges] + [edge[1] for edge in edges]))

				file.write("{}\t{}".format("Dest\\Source", "\t".join(nodes)))

				for dest in nodes:

					file.write("\n{}".format(dest))

					for source in nodes:

						file.write("\t{}".format(edge_counts.get((source, dest), 0)))



	def bootstrap(self, bootstrap_count, sample_size=None):

		for i in range(0, bootstrap_count):

			print("Running bootstrap iteration {}/{}...".format(i + 1, bootstrap_count))

			resampled_aln = self.resample_aln(sample_size)

			(result, modified_inputs) = self.main_analysis(resampled_aln, self.ro_labels)

			self.bs_results_list.append(result[-1])

		#return self.bs_results_list



	def resample_aln(self, sample_size=None):

		if sample_size is None:

			sample_size = len(self.alignment)

		seq_record_list = []

		for i in range(0, sample_size):

			seq_record_list.append(copy.deepcopy(self.ro_aln[random.randint(0, len(self.alignment) - 1)]))

		return MultipleSeqAlignment(seq_record_list)



	def analyze_bootstrap_edges(self):

		main_result = self.main_results[-1]

		bootstrap_results = self.bs_results_list

		bs_edge_counts = {}

		bs_count = len(bootstrap_results)

		for graph in bootstrap_results:

			node_label_dict = {node.obj_dict['name']: node for node in graph.get_nodes()}

			for key in node_label_dict.keys():

				if key == "root":

					node_label_dict[key] = "root"

				else:

					node_label_dict[key] = node_label_dict[key].obj_dict['attributes']['label'].split('\n')[0]

			edges = [(node_label_dict[str(edge.get_source())], node_label_dict[str(edge.get_destination())]) for edge in

					 graph.get_edges()

					 if "color" not in edge.obj_dict['attributes'].keys() or edge.obj_dict['attributes'][

						 'color'] == "green"]

			for edge in edges:

				bs_edge_counts[edge] = bs_edge_counts.get(edge, 0) + 1

		node_label_dict = {node.obj_dict['name']: node for node in main_result.get_nodes()}

		for key in node_label_dict.keys():

			if key == "root":

				node_label_dict[key] = "root"

			else:

				node_label_dict[key] = node_label_dict[key].obj_dict['attributes']['label'].split('\n')[0]

		annotated_edges = []

		for edge in main_result.get_edges():

			if 'color' in edge.obj_dict['attributes'].keys() and edge.obj_dict['attributes']['color'] != "green":

				continue

			source = str(edge.get_source())

			dest = str(edge.get_destination())

			count = bs_edge_counts.get((node_label_dict[source], node_label_dict[dest]), 0)

			if source == "root":

				edge.obj_dict['attributes']['label'] = "({}/{})".format(count, bs_count)

			else:

				edge.obj_dict['attributes']['label'] = edge.obj_dict['attributes']['label'] + "\n({}/{})".format(count,

																												 bs_count)

			annotated_edges.append((node_label_dict[source], node_label_dict[dest]))

		# Add edges displaying alternate bootstrap path weights

		bs_display_thresh = bs_count / 10.0

		node_label_dict = {node_label_dict[key]: key for key in node_label_dict.keys()}

		for edge in bs_edge_counts.keys():

			if edge in annotated_edges or bs_edge_counts[edge] < bs_display_thresh:

				continue

			main_result.add_edge(pydot.Edge(node_label_dict[edge[0]], node_label_dict[edge[1]], color="grey",

											label="{}/{}".format(bs_edge_counts[edge], bs_count), constraint=False))

		return bs_edge_counts



	def tree_to_pydot(self, root, input_aln, labels, comatrix, comatrix_counts):

		graph = pydot.Dot(graph_type='digraph')

		rank = 0

		graph.add_node(pydot.Node(name=root[0], rank=min, label="root"))

		self.pydot_add_children(graph, root, rank + 1, [], input_aln, labels, comatrix, comatrix_counts)

		return graph



	def pydot_add_children(self, graph, node, rank, path, input_aln, labels, comatrix, comatrix_counts):

		child_count = 0

		for child in node[1]:

			if int(round(comatrix[child[0]][child[0]] * self.species_count, 0)) == 0:

				print("Skipping {}...".format(child[0]))

				continue

			new_path = path + [child[0]]

			temp_aln = [seq for seq in input_aln if seq.seq[labels.index(child[0])] == presence_char]

			node_label = "{}\ncount:{}".format(child[0], int(round(comatrix_counts[child[0]][child[0]], 0)))

			if verbose_stats:

				exact_seq = [presence_char if labels[i] in new_path else absence_char for i in range(0, len(labels))]

				exact_seq = "".join(exact_seq)

				hap_exact_count = len([seq for seq in temp_aln if seq.seq[:len(labels)] == exact_seq])

				node_label = node_label + "\nhap_all_count:{}\nhap_exact_count:{}".format(len(temp_aln),

																						  hap_exact_count)

			graph.add_node(

				pydot.Node(name=labels.index(child[0]), rank=rank, label=node_label, path=";".join(new_path)))

			child_count += len(temp_aln)

			if node[0] == "root":

				graph.add_edge(pydot.Edge(node[0], labels.index(child[0])))

			else:

				parent_child_co_occurences = comatrix[child[0]][node[0]]

				graph.add_edge(pydot.Edge(labels.index(node[0]), labels.index(child[0]),

										  label="  {}%".format(round(parent_child_co_occurences*100,1)), penwidth=0.5 + (

								10.0 * (float(comatrix_counts[child[0]][node[0]]) / float(self.species_count)))))

			self.pydot_add_children(graph, child, rank + 1, new_path, temp_aln, labels, comatrix, comatrix_counts)

		if verbose_stats and node[0] != "root":

			pydot_node = graph.get_node(str(labels.index(node[0])))[0]

			pydot_node.obj_dict['attributes']['label'] = pydot_node.obj_dict['attributes'][

															 'label'] + "\nchild_count:{}".format(child_count)



	def highlight_missing(self, graph, comatrix, comatrix_counts, root, labels):

		nodes = graph.get_nodes()

		node_dict = get_node_dict(root)

		link_candidates = set()

		for node in nodes:

			# print(node.obj_dict.keys())

			if node.obj_dict['name'] == "root":

				continue

			node_stats = parse_node_label(node.obj_dict['attributes']['label'])

			if node_stats["count"] > self.species_count * 0.005 and False:

				if node_stats["hap_exact_count"] + node_stats["child_count"] < co_occurence_thresh * node_stats[

					"hap_all_count"]:

					node.obj_dict['attributes']['color'] = "red"

				if node_stats["hap_all_count"] < co_occurence_thresh * node_stats["count"]:

					node.obj_dict['attributes']['color'] = "blue"

			# print(node.obj_dict['attributes'].keys())

			co_freqs = []

			anc_list = node.obj_dict['attributes']['path'].split(";")

			desc_list = get_descendants(node_dict[labels[int(node.obj_dict['name'])]])

			path_list = anc_list + desc_list

			dropped_list = []

			self_freq = comatrix[labels[int(node.obj_dict['name'])]][labels[int(node.obj_dict['name'])]]

			for node_name in [x.obj_dict['name'] for x in nodes]:

				if node_name == "root" or labels[int(node_name)] in desc_list:

					continue

				if labels[int(node_name)] in anc_list:

					if comatrix[labels[int(node.obj_dict['name'])]][labels[int(node_name)]] < 0.5 * self_freq:

						dropped_list.append(labels[int(node_name)])

					continue

				freq = comatrix_counts[labels[int(node.obj_dict['name'])]][labels[int(node_name)]]

				co_freqs.append((node_name, freq))

			co_freqs = [(labels[int(x[0])], x[1]) for x in co_freqs]

			co_freqs.sort(key=lambda tup: tup[1], reverse=True)

			co_freq_top_count = 3

			if "color" in node.obj_dict['attributes'].keys():

				node.obj_dict['attributes']['penwidth'] = 5.0

				co_freq_top_count = 3

			for i in range(0, min(co_freq_top_count, len(co_freqs))):

				node.obj_dict['attributes']['label'] = node.obj_dict['attributes']['label'] + "\n{}: {}".format(

					co_freqs[i][0], co_freqs[i][1])

			if len(dropped_list) > 0:

				node.obj_dict['attributes']['label'] = node.obj_dict['attributes']['label'] + "\n{}: {}".format(

					"Dropped", ",".join(dropped_list))

			missing_count_threshold = self.species_count * 0.01

			# if co_freqs[0][1] > min(missing_count_threshold, node_stats["count"] * 0.25):

			# if len(co_freqs) > 0 and (co_freqs[0][1] > missing_count_threshold or "color" in node.obj_dict['attributes'].keys()):

			if len(co_freqs) > 0 and co_freqs[0][1] > missing_count_threshold:

				i = 0

				linked_node = co_freqs[0][0]

				for co_freq in co_freqs:

					if co_freq[1] > 0.8 * co_freqs[0][1]:

						i += 1

				co_freqs = co_freqs[:i]

				for i in range(0, len(co_freqs)):

					for co_freq in co_freqs:

						if co_freq[0] in [x[0] for x in node_dict[linked_node][1]]:

							linked_node = co_freq[0]

				link_candidates.add((labels[int(node.obj_dict['name'])], linked_node))

		link_candidates = list(link_candidates)

		links = []

		for link_candidate in link_candidates:

			children_L = [x[0] for x in node_dict[link_candidate[0]][1]]

			children_R = [x[0] for x in node_dict[link_candidate[1]][1]]

			child_links = [(link_candidate[0], x) for x in children_R] + [(link_candidate[1], x) for x in children_L]

			child_links = child_links + [(x, link_candidate[0]) for x in children_R] + [(x, link_candidate[1]) for x in

																						children_L]

			has_child_link = False

			for child_link in child_links:

				if child_link in link_candidates:

					has_child_link = True

					break

			if not has_child_link:

				links.append(link_candidate)

		for link in links:

			# graph.add_edge(pydot.Edge(labels.index(link[0]), labels.index(link[1]), penwidth=3.0, color="#ff00005f", constraint=False))

			graph.add_edge(pydot.Edge(labels.index(link[1]), labels.index(link[0]), penwidth=3.0, color="#ff00005f",

									  constraint=False))



	def shade_by_io_diff(self, input_graph):

		# color_scale = list(Color("red").range_to(Color("black"), 50)) + list(Color("black").range_to(Color("blue"), 50))

		# line_div_scale = 500.0

		line_div_scale = 0.03 * self.species_count

		# color_div_scale = 50.0

		color_div_scale = 0.003 * self.species_count

		# diff_threshold = 500.0

		diff_threshold = 0.03 * self.species_count

		for node in input_graph.get_nodes():

			if node.obj_dict['name'] == "root":

				continue

			node_data = parse_node_label(node.obj_dict['attributes']['label'])

			diff = node_data["count"] - (node_data["child_count"] + node_data["hap_exact_count"])

			if abs(diff) < diff_threshold:

				continue

			linewidth = 1.0 + min(9.0, (abs(diff) / line_div_scale))

			color_coord = int(min(max(diff / color_div_scale, -50.0), 49.0) + 50.0)

			node.obj_dict["attributes"]["penwidth"] = linewidth

			# node.obj_dict["attributes"]["color"] = str(color_scale[color_coord])

			if diff < 0:

				node.obj_dict["attributes"]["color"] = "red"

			elif diff > 0:

				node.obj_dict["attributes"]["color"] = "blue"



	def reverse_overpopulated(self, input_graph, initial_aln, labels):

		nodes = input_graph.get_nodes()

		new_aln = copy.deepcopy(initial_aln)

		rev_pos_list = []

		edges = [(edge.get_source(), edge.get_destination()) for edge in input_graph.get_edges()]

		target = ""

		max_dropped = 0

		for node in nodes:

			if node.obj_dict['name'] == "root":

				continue

			label_string = node.obj_dict['attributes']['label']

			node_data = parse_node_label(label_string)

			#print(node.obj_dict['name'])

			# if "17801C" in labels[int(node.obj_dict['name'])]:

			# 	print("Node:{}\nSpecies count:{}\nChild count:{}\nCount:{}\n".format(labels[int(node.obj_dict['name'])], self.species_count, node_data["child_count"], node_data["count"]))

			if node_data["child_count"] > node_data["count"] + (self.species_count * 0.01):

				children = [edge[1] for edge in edges if edge[0] == int(node.obj_dict['name'])]

				# target = ""

				# max_dropped = 0

				for child in children:

					child = input_graph.get_node(str(child))[0]

					if labels[int(child.obj_dict['name'])] in self.comped_nodes:

						continue

					child_data = parse_node_label(child.obj_dict['attributes']['label'])

					dropped = child_data["hap_all_count"] - (child_data["child_count"] + child_data["hap_exact_count"])

					if dropped < child_data["hap_all_count"] * 0.5 or child_data["count"] < self.species_count * 0.1:

						continue

					if dropped > max_dropped:

						max_dropped = dropped

						target = child.obj_dict['name']

				if max_dropped > 0:

					pass

					print("Binary complementing {}...".format(labels[int(target)]))

					# print("Binary complementing {}...".format(label_string.split("\n")[0].strip()))

					# print("Node:{}\nSpecies count:{}\nChild count:{}\nCount:{}\n".format(labels[int(node.obj_dict['name'])], self.species_count, node_data["child_count"], node_data["count"]))

					# rev_pos_list.append(int(target))

			# new_aln[int(node.obj_dict['name'])].seq = get_binary_complement(new_aln[int(node.obj_dict['name'])].seq)

		if target != "":

			rev_pos_list.append(int(target))

			binary_complement_positions(new_aln, rev_pos_list)

			# print([labels[int(x)] for x in self.reversed_set])

		return new_aln, rev_pos_list



	def spr_pass(self, graph, comatrix, comatrix_counts, root, labels):

		nodes = graph.get_nodes()

		node_dict = get_node_dict(root)

		anc_lists = {}

		desc_lists = {}

		path_lists = {}

		co_lists = {}

		depths = {}

		regrafted = []

		new_edges = []

		nodes_by_depth = {}

		sorted_node_list = []

		# print("labels: {}".format(labels))

		# print("init_nodes: {}".format(self.initial_nodes))

		for node in nodes:

			# print(node.obj_dict.keys())

			node_name = node.obj_dict['name']

			if node_name == "root":

				continue

			node_data = parse_node_label(node.obj_dict['attributes']['label'])

			anc_lists[node_name] = node.obj_dict['attributes']['path'].split(";")

			desc_lists[node_name] = get_descendants(node_dict[labels[int(node.obj_dict['name'])]])

			path_lists[node_name] = anc_lists[node_name] + desc_lists[node_name]

			co_lists[node_name] = [key for key in comatrix[labels[int(node.obj_dict['name'])]].keys() if

								   comatrix[labels[int(node.obj_dict['name'])]][key] >= co_occurence_thresh]

			# co_lists[node_name] = [key for key in comatrix[labels[int(node.obj_dict['name'])]].keys() if

			# 					   comatrix[labels[int(node.obj_dict['name'])]][key] >= 0.51]

			# print(node_name)

			# print(labels[int(node_name)])

			# print(co_lists[node_name])

			depths[node_name] = get_depth(root, labels[int(node.obj_dict['name'])])

			nodes_by_depth[depths[node_name]] = nodes_by_depth.get(depths[node_name], []) + [node]

		# for key in sorted(nodes_by_depth.keys()):

		# 	sorted_node_list += nodes_by_depth[key]

		for node in nodes:

			if node.obj_dict['name'] == "root":

				continue

			sorted_node_list.append((node, comatrix_counts[labels[int(node.obj_dict['name'])]][labels[int(node.obj_dict['name'])]]))

		sorted_node_list.sort(key=lambda tup: tup[1], reverse=True)

		sorted_node_list = [x[0] for x in sorted_node_list]

		for node in sorted_node_list:

			node_name = node.obj_dict['name']

			if node_name == "root" or labels[int(node_name)] in regrafted or (self.params.lock_init_tree and labels[int(node_name)] in self.initial_nodes):

				continue

			desc_list = get_descendants(node_dict[labels[int(node.obj_dict['name'])]])

			# co_path_overlap = sum([1 for mut in co_lists[node_name] if mut in path_lists[node_name]])

			co_path_overlap = sum([1 for mut in co_lists[node_name] if mut in anc_lists[node_name]])

			# if float(co_path_overlap)/float(len(co_lists[node_name])) > 0.5:

			# if float(co_path_overlap) / float(len(co_lists[node_name])) > 0.97:

			if co_path_overlap == len(co_lists[node_name]):

				continue

			new_parent = ""

			max_match_score = 0

			max_depth = 0

			for target_node_name in anc_lists.keys():

				if labels[int(target_node_name)] in desc_list:

					continue

				co_anc = sum([1 for mut in co_lists[node_name] if mut in anc_lists[target_node_name]])

				mismatch_anc = sum([1 for mut in anc_lists[target_node_name] if mut not in co_lists[node_name]])

				match_score = co_anc - (mismatch_anc * self.params.spr_mismatch_weight)

				# if node_name == '45':

				# 	print("{} - {}: {}".format(labels[int(node_name)], labels[int(target_node_name)], match_score))

				if match_score > max_match_score:

					max_match_score = match_score

					max_depth = depths[target_node_name]

					new_parent = target_node_name

				elif match_score == max_match_score and depths[target_node_name] > max_depth:

					max_depth = depths[target_node_name]

					new_parent = target_node_name

			if vlevel > 2:

				print("Already regrafted: {}".format(",".join(regrafted)))

			if new_parent not in ["", node_name]:

				if vlevel >= 2:

					print("Regrafting {} as child of {}...".format(labels[int(node_name)], labels[int(new_parent)]))

				regrafted.extend(desc_list)

				regraft(graph, root, labels, labels[int(node_name)], labels[int(new_parent)])

				new_edges.append((labels[int(new_parent)], labels[int(node_name)]))

		return new_edges



	def add_timing_labels(self, graph, var_days):

		nodes = graph.get_nodes()

		if self.params.all_times is not None:

			with open(self.params.all_times, 'r') as times_file:

				all_times = sorted([int(line.strip()) for line in times_file])

		for node in nodes:

			if node.obj_dict['name'] == "root":

				continue

			label = node.obj_dict['attributes']['label']

			label_lines = label.split('\n')

			name = label_lines[0].strip()

			if "_comp" in name:

				name = name.replace("_comp", "")

				if self.params.all_times is None:

					obs_list = var_days[name]

				else:

					obs_list = copy.deepcopy(all_times)

					for time in var_days[name]:

						obs_list.remove(time)

			else:

				obs_list = var_days[name]

			obs_cnt = len(obs_list) - 1

			quarts = [str(obs_list[int(x * 0.25 * obs_cnt)]) for x in range(0, 5)]

			# new_label_lines = [label_lines[0], "seen_qrts: {}".format(",".join(quarts))]

			new_label_lines = [label_lines[0], "first_seen: {}".format(quarts[0])]

			new_label_lines = new_label_lines + label_lines[1:]

			node.obj_dict['attributes']['label'] = "\n".join(new_label_lines)



	def assign_seqs_to_nodes(self, input_aln, labels, graph):

		# Get set of ancestral mutations for each node

		ancestors = graph_to_ancestors_dict(graph)

		# Make sequence for each node

		node_list = []

		node_seq_list = []

		# print(", ".join(labels))

		# print(", ".join(ancestors.keys()))

		for key in ancestors.keys():

			node_list.append(key)

			node_seq_list.append(numpy.array([1 if label in ancestors[key] else 0 for label in labels]))

		# Do matrix math thing to identify node that maximizes mutation matches and then minimizes mutation mismatches for each sequence

		zeros = numpy.zeros(len(labels))

		binary_seq_list = [numpy.array([-1 if sequence.seq[i] == presence_char else 0 for i in range(0, len(labels))])

						   for sequence in input_aln]

		binary_seq_set = list(numpy.unique(binary_seq_list, axis=0))

		binary_seq_set_strings = ["".join([str(z) for z in y]) for y in binary_seq_set]

		binary_seq_best_nodes = []

		binary_seq_best_counts = []

		for binary_seq in binary_seq_set:

			diff_counts = []  # List of (missed_mut_count, extra_mut_count, matched_mut_count) tuples

			for node_seq in node_seq_list:

				row = numpy.add(node_seq, binary_seq)

				missed_mut_count = abs(sum(numpy.minimum(row, zeros)))

				extra_mut_count = abs(sum(numpy.maximum(row, zeros)))

				matched_mut_count = abs(sum(binary_seq)) - missed_mut_count

				diff_counts.append((missed_mut_count, extra_mut_count, matched_mut_count))

			min_missed = min([x[0] for x in diff_counts])

			min_extra = min([x[1] for x in diff_counts if x[0] == min_missed])

			max_matched = max([x[2] for x in diff_counts if (x[0], x[1]) == (min_missed, min_extra)])

			if diff_counts.count((min_missed, min_extra, max_matched)) != 1:

				best_node = []

				if vlevel > 3:

					print("Found sequence with best fit at {} nodes...".format(

						diff_counts.count((min_missed, min_extra, max_matched))))

					print(binary_seq)

				for i in range(0, len(diff_counts)):

					if diff_counts[i] == (min_missed, min_extra, max_matched):

						best_node.append(node_list[i])

						if vlevel > 3:

							print(node_list[i])

							print(node_seq_list[i])

			else:

				best_node = [node_list[diff_counts.index((min_missed, min_extra, max_matched))]]

			binary_seq_best_nodes.append((best_node, (min_missed, min_extra, max_matched)))

		best_fit_nodes = []

		for i in range(0, len(input_aln)):

			# print(binary_seq_set_strings.index("".join([str(z) for z in binary_seq_list[i]])))

			# print(binary_seq_best_nodes[binary_seq_set_strings.index("".join([str(z) for z in binary_seq_list[i]]))])

			# print(binary_seq_set.index(binary_seq_list[i]))

			best_fit_nodes.append((input_aln[i].id, ",".join([str(abs(z)) for z in binary_seq_list[i]]),

								   binary_seq_best_nodes[

									   binary_seq_set_strings.index("".join([str(z) for z in binary_seq_list[i]]))][0],

								   binary_seq_best_nodes[

									   binary_seq_set_strings.index("".join([str(z) for z in binary_seq_list[i]]))][1]))

		print("Totals -- Matches:{}, Misses:{}, Extras:{}".format(sum([x[3][2] for x in best_fit_nodes]), sum([x[3][0] for x in best_fit_nodes]), sum([x[3][1] for x in best_fit_nodes])))

		if self.params.stats_file is not None:

			with open(self.params.stats_file,'a+') as file:

				file.write("{}\t{}\t{}\t{}\n".format(sum([x[3][2] for x in best_fit_nodes]), sum([x[3][0] for x in best_fit_nodes]), sum([x[3][1] for x in best_fit_nodes]), ",".join(list(self.comped_nodes))))

		with open(os.path.join(self.params.output, "seq_assign_Test.txt"), 'w') as file:

			file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("seq_id", "best_node", "miss_cnt", "extra_cnt", "match_cnt",

														 "binary_seq({})".format(",".join(labels))))

			for val in best_fit_nodes:

				file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(val[0], ",".join(val[2]), val[3][0], val[3][1], val[3][2],

															 val[1]))



	def detect_reversals(self, input_aln, input_labels, matrix_cache_basename):

		aln = copy.deepcopy(input_aln)

		labels = copy.deepcopy(input_labels)

		print("Calculating comatrix...")

		(comatrix, comatrix_counts) = build_comatrix(aln, labels)

		sorted_labels = []

		for label in labels:

			sorted_labels.append((label, comatrix_counts[label][label]))

		sorted_labels.sort(key=lambda tup: tup[1], reverse=True)

		for i in range(0, len(aln)):

			seq_str = aln[i]

			new_seq = Seq("".join([seq_str[labels.index(label[0])] for label in sorted_labels]), aln[i].seq.alphabet)

			aln[i].seq = new_seq

		labels = [x[0] for x in sorted_labels]

		votes = []

		print("Detecting reversals...")

		with open("{}_rev.txt".format(matrix_cache_basename), 'w') as outfile2:

			with open("{}_fwd.txt".format(matrix_cache_basename), 'w') as outfile:

				outfile.write("\t{}".format("\t".join(labels)))

				outfile2.write("\t{}".format("\t".join(labels)))

				for i in range(0, len(labels)):

					flip_votes = 0

					label2 = labels[i]

					outfile.write("\n{}".format(label2))

					outfile2.write("\n{}".format(label2))

					for j in range(0, len(labels)):

						label1 = labels[j]

						both = comatrix_counts[label1][label2]

						only1 = comatrix_counts[label1][label1] - both

						only2 = comatrix_counts[label2][label2] - both

						neither = self.species_count - (only1 + only2 + both)

						rev2_both = only1

						rev2_neither = only2

						rev2_only2 = neither

						rev2_only1 = both

						min_fwd = min(only1, both, only2) / float(both + only2)

						min_rev = min(rev2_only1, rev2_both, rev2_only2) / float(rev2_both + rev2_only2)

						outfile.write("\t{}".format(min(only1, both, only2)))

						outfile2.write("\t{}".format(min(rev2_only1, rev2_both, rev2_only2)))

						# if j < 3:

						# 	print("{}/{}:({}, {}, {})\t{}/{}_rev:({}, {}, {})".format(label1, label2, only1, both, only2, label1, label2, rev2_only1, rev2_both, rev2_only2))

						# if (min_rev < min_fwd and label1 not in [val[0] for val in votes if val[1] > val[2]*0.5]) or (min_rev > min_fwd and label1 in [val[0] for val in votes if val[1] > val[2]*0.5]):

						if min_rev < min_fwd:

							flip_votes += 1

					votes.append((label2, flip_votes, len(labels)))

					if flip_votes > 0.5*len(labels):

						print("Voted {}/{} times to flip {}...".format(flip_votes, len(labels), label2))

				for val in votes:

					print("{}:{}/{}".format(val[0], val[1], val[2]))

		max_votes = max([val[1] for val in votes])

		print("Complementing: {}".format([val[0] for val in votes if val[1] == max_votes]))

		binary_complement_positions(aln, [labels.index(val[0]) for val in votes if val[1] == max_votes])

		return aln, labels





	def detect_reversals_new(self, input_aln, input_labels, matrix_cache_basename):

		aln = copy.deepcopy(input_aln)

		labels = copy.deepcopy(input_labels)

		print("Calculating comatrix...")

		(comatrix, comatrix_counts) = build_comatrix(aln, labels)

		sorted_labels = []

		for label in labels:

			sorted_labels.append((label, comatrix_counts[label][label]))

		sorted_labels.sort(key=lambda tup: tup[1], reverse=True)

		for i in range(0, len(aln)):

			seq_str = aln[i]

			new_seq = Seq("".join([seq_str[labels.index(label[0])] for label in sorted_labels])) #, aln[i].seq.alphabet)

			aln[i].seq = new_seq

		labels = [x[0] for x in sorted_labels]

		votes = []

		scores = []

		min_fwd_cnt_dict = {}

		min_rev_cnt_dict = {}

		print("Detecting reversals...")

		with open("{}_rev.txt".format(matrix_cache_basename), 'w') as outfile2:

			with open("{}_fwd.txt".format(matrix_cache_basename), 'w') as outfile:

				outfile.write("\t{}".format("\t".join(labels)))

				outfile2.write("\t{}".format("\t".join(labels)))

				min_total = 0

				for i in range(0, len(labels)):

					flip_votes = 0

					flip_score = 0

					min_fwd_cnt = 0

					total_fwd_cnt = 0

					min_rev_cnt = 0

					total_rev_cnt = 0

					label2 = labels[i]

					min_fwd_cnt_dict[label2] = {}

					min_rev_cnt_dict[label2] = {}

					outfile.write("\n{}".format(label2))

					outfile2.write("\n{}".format(label2))

					for j in range(0, len(labels)):

						label1 = labels[j]

						both = comatrix_counts[label1][label2]

						only1 = comatrix_counts[label1][label1] - both

						only2 = comatrix_counts[label2][label2] - both

						neither = self.species_count - (only1 + only2 + both)

						rev2_both = only1

						rev2_neither = only2

						rev2_only2 = neither

						rev2_only1 = both

						# min_fwd = min(only1, both, only2) / float(both + only2)

						# min_rev = min(rev2_only1, rev2_both, rev2_only2) / float(rev2_both + rev2_only2)

						min_fwd = min(only1, both, only2) / comatrix_counts[label2][label2]

						min_rev = min(rev2_only1, rev2_both, rev2_only2) / (self.species_count - comatrix_counts[label2][label2])

						outfile.write("\t{}".format(min(only1, both, only2)))

						outfile2.write("\t{}".format(min(rev2_only1, rev2_both, rev2_only2)))

						# if j < 3:

						# 	print("{}/{}:({}, {}, {})\t{}/{}_rev:({}, {}, {})".format(label1, label2, only1, both, only2, label1, label2, rev2_only1, rev2_both, rev2_only2))

						# if (min_rev < min_fwd and label1 not in [val[0] for val in votes if val[1] > val[2]*0.5]) or (min_rev > min_fwd and label1 in [val[0] for val in votes if val[1] > val[2]*0.5]):

						flip_score += min_fwd - min_rev

						min_fwd_cnt += min(only1, both, only2)

						min_fwd_cnt_dict[label2][label1] = min(only1, both, only2)

						min_rev_cnt += min(rev2_only1, rev2_both, rev2_only2)

						min_rev_cnt_dict[label2][label1] = min(rev2_only1, rev2_both, rev2_only2)

						total_fwd_cnt += float(both + only2)

						total_rev_cnt += float(rev2_both + rev2_only2)

						# total_fwd_cnt += float(only1 + both + only2)

						# total_rev_cnt += float(rev2_only1 + rev2_both + rev2_only2)

						if min_rev < min_fwd * self.params.flip_vote_thresh:

							flip_votes += 1

					votes.append((label2, flip_votes, len(labels)))

					# flip_score = (min_fwd_cnt/total_fwd_cnt) - (min_rev_cnt/total_rev_cnt)

					# flip_score = (min_fwd_cnt - min_rev_cnt)/min_fwd_cnt

					flip_score = (min_fwd_cnt/comatrix_counts[label2][label2]) - (min_rev_cnt / (self.species_count - comatrix_counts[label2][label2]))

					min_total += min_fwd_cnt

					scores.append((label2, flip_score, min_fwd_cnt, min_rev_cnt))

					# if flip_votes > 0.5*len(labels):

					# 	print("Voted {}/{} times to flip {} with a score of {}".format(flip_votes, len(labels), label2, round(flip_score,4)))

					if flip_score > 0 and flip_votes > len(labels) * 0.5:

						print("{}/{} votes to flip {}, score of {} ({}/{} -> {}/{})".format(flip_votes, len(labels), label2, round(flip_score, 4),

																										min_fwd_cnt, comatrix_counts[label2][label2],

																										min_rev_cnt, self.species_count - comatrix_counts[label2][label2]))

				# for val in votes:

				# 	print("{}:{}/{}".format(val[0], val[1], val[2]))

				print("Total flip votes: {}".format(sum([val[1] for val in votes])))

				print("Total fwd min count: {}".format(min_total))

				print("Total comatrix counts: {}".format(sum([sum(val.values()) for val in comatrix_counts.values()])))

				pdx_ratio = min_total / sum([sum(val.values()) for val in comatrix_counts.values()])

				print("Paradox ratio: {}".format(round(pdx_ratio, 5)))



		# max_votes = max([val[1] for val in votes])

		# print("Complementing: {}".format([val[0] for val in votes if val[1] == max_votes]))

		# for val in votes:

		# 	if val[1] == max_votes:

		# 		print("Complementing: {} with score {}...".format(val[0], val[1]))

		# binary_complement_positions(aln, [labels.index(val[0]) for val in votes if val[1] == max_votes])

		min_fwd_cnt_list = []

		for key1 in min_fwd_cnt_dict.keys():

			min_fwd_cnt_list.extend(min_fwd_cnt_dict[key1].values())

		min_fwd_cnt_list.sort(reverse=True)

		keys = list(min_fwd_cnt_dict.keys())

		for i in range(0, len(keys)):

			key1 = keys[i]

			for j in range(i, len(keys)):

				key2 = keys[j]

				if min_fwd_cnt_dict[key1][key2] in min_fwd_cnt_list[0:5]:

					print("High paradox count({}) for labels {} and {}".format(min_fwd_cnt_dict[key1][key2], key1, key2))

					both = comatrix_counts[key1][key2]

					only1 = comatrix_counts[key1][key1] - both

					only2 = comatrix_counts[key2][key2] - both

					print("{}:{}\t{}:{}\tBoth:{}".format(key1, only1, key2, only2, both))

					for val in scores:

						if val[0] in [key1, key2]:

							print("{} flip score: {}".format(val[0], val[1]))

		comp_list = []

		# max_score = max([val[1] for val in scores])

		# if max_score > 0:

		# 	print("### Complementing: {} ###".format([val[0] for val in scores if val[1] == max_score]))

		# 	for val in scores:

		# 		if val[1] == max_score:

		# 			print("Complementing: {} with score {} ({} - {})...".format(val[0], round(val[1], 4), val[2], val[3]))

		# 			print("Old count: {}\tNew count: {}\tReduction: {}".format(comatrix_counts[val[0]][val[0]], self.species_count - comatrix_counts[val[0]][val[0]], round(2 - (self.species_count/comatrix_counts[val[0]][val[0]]), 4)))

		# 			mp_list = sorted(list(min_fwd_cnt_dict[val[0]].values()), reverse=True)[0:5]

		# 			mp_keys = [key for key in min_fwd_cnt_dict[val[0]].keys() if min_fwd_cnt_dict[val[0]][key] in mp_list]

		# 			for key in mp_keys:

		# 				print("{} has high paradox count {} with {}".format(val[0], min_fwd_cnt_dict[val[0]][key], key))

		# 			#print("{} has max paradox count {} with {}".format(val[0], mp_cnt, mp_keys))

		# 	binary_complement_positions(aln, [labels.index(val[0]) for val in scores if val[1] == max_score])

		max_votes = max([val[1] for val in votes if sum(min_fwd_cnt_dict[val[0]].values())/(comatrix_counts[val[0]][val[0]]) > sum(min_rev_cnt_dict[val[0]].values())/(self.species_count - comatrix_counts[val[0]][val[0]])])

		if max_votes > len(input_labels) * self.params.flip_pass_thresh:

			print("### Potentially complementing: {} with {} votes.###".format([val[0] for val in votes if val[1] == max_votes], max_votes))

			for val in votes:

				if val[1] == max_votes and sum(min_fwd_cnt_dict[val[0]].values())/(comatrix_counts[val[0]][val[0]]) > sum(min_rev_cnt_dict[val[0]].values())/(self.species_count - comatrix_counts[val[0]][val[0]]):

					comp_list.append(val[0])

					print("Complementing: {} with {} votes ({} - {})...".format(val[0], round(val[1], 4), sum(min_fwd_cnt_dict[val[0]].values()), sum(min_rev_cnt_dict[val[0]].values())))

					print("Old count: {}\tNew count: {}".format(comatrix_counts[val[0]][val[0]], self.species_count - comatrix_counts[val[0]][val[0]]))

					mp_list = sorted(list(min_fwd_cnt_dict[val[0]].values()), reverse=True)[0:5]

					mp_keys = [key for key in min_fwd_cnt_dict[val[0]].keys() if min_fwd_cnt_dict[val[0]][key] in mp_list]

					for key in mp_keys:

						print("{} has high paradox count {} with {}".format(val[0], min_fwd_cnt_dict[val[0]][key], key))

					#print("{} has max paradox count {} with {}".format(val[0], mp_cnt, mp_keys))

			binary_complement_positions(aln, [labels.index(val) for val in comp_list])

		else:

			print("Found no candidates for complementing.")

		return aln, labels, comp_list





	def main_analysis(self, input_aln=None, labels=None, var_days=None):

		start_time = datetime.datetime.now()

		if input_aln is None:

			input_aln = copy.deepcopy(self.alignment)

		if labels is None:

			labels = copy.deepcopy(self.mut_labels)

		if var_days is None:

			var_days = copy.deepcopy(self.mut_days)

		if cache_bootstrap_aln:

			aln_cache_filename = "aln_cache_test.fas"

			labels_cache_filename = "label_cache_test.txt"

			with open(os.path.join(self.params.output, aln_cache_filename), 'w') as aln_cache_file:

				AlignIO.write(input_aln, aln_cache_file, "fasta")

			with open(os.path.join(self.params.output, labels_cache_filename), 'w') as labels_cache_file:

				labels_cache_file.write("label\tdays\n")

				for label in labels:

					labels_cache_file.write("{}\t{}\n".format(label, ",".join([str(day) for day in var_days[label]])))

		results_list = []

		temp_aln = copy.deepcopy(input_aln)

		temp_labels = copy.deepcopy(labels)

		comp_list = []

		i = 0

		while True:

			i += 1

			(temp_aln, temp_labels, temp_comp_list) = self.detect_reversals_new(temp_aln, temp_labels, os.path.join(self.params.output, "rev_detection_test_step{}".format(i)))

			if len(temp_comp_list) == 0:

				if i > 1:

					input_aln = temp_aln

					labels = temp_labels

				break

			comp_list.extend(temp_comp_list)

		# exit()

		self.comped_nodes.update(comp_list)

		(comatrix, comatrix_counts) = build_comatrix(input_aln, labels)

		sorted_labels = []

		for label in labels:

			sorted_labels.append((label, comatrix_counts[label][label]))

		sorted_labels.sort(key=lambda tup: tup[1], reverse=True)

		for i in range(0, len(input_aln)):

			seq_str = input_aln[i]

			new_seq = Seq("".join([seq_str[labels.index(label[0])] for label in sorted_labels]))#, input_aln[i].seq.alphabet)

			input_aln[i].seq = new_seq

		labels = [x[0] for x in sorted_labels]

		if self.params.initial_tree is None:

			root_node = build_tree(labels, comatrix)

		else:

			root_node = update_tree(labels, comatrix, self.initial_tree)

		graph = self.tree_to_pydot(root_node, input_aln, labels, comatrix, comatrix_counts)

		if not self.is_bootstrapping:

			graph.write_png(os.path.join(self.params.output, output_basename + "_{}.png".format("TESTING")))

			print("Drawing initial tree...")

		self.highlight_missing(graph, comatrix, comatrix_counts, root_node, labels)

		reset_highlighting(graph)

		self.shade_by_io_diff(graph)

		results_list.append(graph)

		# Excess children complementing pass

		# Preset ro values to current values in case no nodes are reversed

		ro_comatrix = comatrix

		ro_comatrix_counts = comatrix_counts

		ro_root_node = root_node

		ro_graph = graph

		ro_aln = input_aln

		if not self.is_bootstrapping:

			#while len(reversed_list) > 0:

			while True:

				if self.params.disable_graph_flipping:

					reversed_list = []

					if self.params.bootstrap_reps > 0:

						self.ro_aln = copy.deepcopy(ro_aln)

						self.ro_labels = copy.deepcopy(labels)

				else:

					(ro_aln, reversed_list) = self.reverse_overpopulated(ro_graph, ro_aln, labels)

				if len(reversed_list) < 1:

					break

				(ro_comatrix, ro_comatrix_counts) = build_comatrix(ro_aln, labels)

				self.comped_nodes.update([labels[x] for x in reversed_list])

				sorted_labels = []

				for label in labels:

					sorted_labels.append((label, ro_comatrix_counts[label][label]))

				sorted_labels.sort(key=lambda tup: tup[1], reverse=True)

				for i in range(0, len(ro_aln)):

					seq_str = ro_aln[i]

					new_seq = Seq("".join([seq_str[labels.index(label[0])] for label in sorted_labels])) #,

								 # ro_aln[i].seq.alphabet)

					ro_aln[i].seq = new_seq

				labels = [x[0] for x in sorted_labels]

				ro_root_node = build_tree(labels, ro_comatrix)

				ro_graph = self.tree_to_pydot(ro_root_node, ro_aln, labels, ro_comatrix, ro_comatrix_counts)

				self.highlight_missing(ro_graph, ro_comatrix, ro_comatrix_counts, ro_root_node, labels)

				reset_highlighting(ro_graph)

				self.shade_by_io_diff(ro_graph)

				for reversed_pos_label in self.comped_nodes:

					reversed_pos = labels.index(reversed_pos_label)

					node = ro_graph.get_node(str(reversed_pos))[0]

					node_data = node.obj_dict['attributes']['label'].split("\n")

					node_data[0] = node_data[0] + "_comp"

					node.obj_dict['attributes']['label'] = "\n".join(node_data)

					# if "color" not in node.obj_dict['attributes'].keys():

					if True:

						node.obj_dict['attributes']['color'] = "green"

						node.obj_dict['attributes']['penwidth'] = 3.0

				results_list.append(ro_graph)

			self.ro_aln = copy.deepcopy(ro_aln)

			self.ro_labels = copy.deepcopy(labels)

		else:

			# input_aln already set to saved ro_aln and ro_aln_labels by call to main_analysis

			ro_aln = input_aln

			# fetch comped_nodes from self attribute

			for reversed_pos_label in self.comped_nodes:

				reversed_pos = labels.index(reversed_pos_label)

				node = ro_graph.get_node(str(reversed_pos))[0]

				node_data = node.obj_dict['attributes']['label'].split("\n")

				node_data[0] = node_data[0] + "_comp"

				node.obj_dict['attributes']['label'] = "\n".join(node_data)

				# if "color" not in node.obj_dict['attributes'].keys():

				if True:

					node.obj_dict['attributes']['color'] = "green"

					node.obj_dict['attributes']['penwidth'] = 3.0

			results_list.append(ro_graph)

		# Add comp annotations to initial result graph

		for node in results_list[0].get_nodes():

			if node.obj_dict["name"] == "root":

				continue

			node_data = node.obj_dict['attributes']['label'].split("\n")

			if node_data[0] in self.comped_nodes:

				node_data[0] = node_data[0] + "_comp"

				node.obj_dict['attributes']['label'] = "\n".join(node_data)

				node.obj_dict['attributes']['color'] = "green"

				node.obj_dict['attributes']['penwidth'] = 3.0

		with open(os.path.join(self.params.output, "COI_matrix.txt"), 'w') as file:

			keys = ro_comatrix.keys()

			file.write("child\\parent\t{}".format("\t".join(keys)))

			for key1 in keys:

				file.write("\n{}".format(key1))

				for key2 in keys:

					file.write("\t{}".format(round(ro_comatrix[key1][key2], 3)))

		with open(os.path.join(self.params.output, "co_counts_matrix.txt"), 'w') as file:

			keys = ro_comatrix.keys()

			file.write("child\\parent\t{}".format("\t".join(keys)))

			for key1 in keys:

				file.write("\n{}".format(key1))

				for key2 in keys:

					file.write("\t{}".format(ro_comatrix_counts[key1][key2], 3))

		spr_graph = self.tree_to_pydot(ro_root_node, ro_aln, labels, ro_comatrix, ro_comatrix_counts)

		grafted_edges = self.spr_pass(spr_graph, ro_comatrix, ro_comatrix_counts, ro_root_node, labels)

		spr_graph = self.tree_to_pydot(ro_root_node, ro_aln, labels, ro_comatrix, ro_comatrix_counts)

		self.highlight_missing(spr_graph, ro_comatrix, ro_comatrix_counts, ro_root_node, labels)

		# reset_highlighting(spr_graph)

		self.shade_by_io_diff(spr_graph)

		for reversed_pos_label in self.comped_nodes:

			reversed_pos = labels.index(reversed_pos_label)

			node = spr_graph.get_node(str(reversed_pos))[0]

			node_data = node.obj_dict['attributes']['label'].split("\n")

			node_data[0] = node_data[0] + "_comp"

			node.obj_dict['attributes']['label'] = "\n".join(node_data)

			if True:

				node.obj_dict['attributes']['color'] = "green"

				node.obj_dict['attributes']['penwidth'] = 3.0

		for edge in grafted_edges:

			edge_obj = spr_graph.get_edge(labels.index(edge[0]), labels.index(edge[1]))[0]

			edge_obj.obj_dict['attributes']['color'] = "green"

			edge_obj.obj_dict['attributes']['penwidth'] = 3.0

		# highlight_missing(spr_graph, ro_comatrix, ro_comatrix_counts,, ro_root_node, labels)

		results_list.append(spr_graph)

		for i in range(0, len(labels)):

			if labels[i] in self.comped_nodes:

				labels[i] = labels[i] + "_comp"

		if vlevel > 1:

			print("Main analysis time:{}".format(datetime.datetime.now() - start_time))

		return results_list, (ro_aln, labels, var_days)





JobBase(args)




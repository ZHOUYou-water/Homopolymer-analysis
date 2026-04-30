import pandas as pd

# 文件路径
names_file = "~/Desktop/names.dmp"
nodes_file = "~/Desktop/nodes.dmp"

# 读取 names.dmp
names = pd.read_csv(
    names_file,
    sep="|",
    header=None,
    usecols=[0, 1, 3],
    names=["tax_id", "name_txt", "name_class"],
    dtype=str
)
names["tax_id"] = names["tax_id"].str.strip()
names["name_txt"] = names["name_txt"].str.strip()
names["name_class"] = names["name_class"].str.strip()

# 只保留 scientific name
scientific_names = names[names["name_class"] == "scientific name"].drop(columns=["name_class"])

# 读取 nodes.dmp
nodes = pd.read_csv(
    nodes_file,
    sep="|",
    header=None,
    usecols=[0, 1, 2],
    names=["tax_id", "parent_tax_id", "rank"],
    dtype=str
)
nodes["tax_id"] = nodes["tax_id"].str.strip()
nodes["parent_tax_id"] = nodes["parent_tax_id"].str.strip()
nodes["rank"] = nodes["rank"].str.strip()

# 合并 rank 和 scientific name
taxonomy = pd.merge(nodes, scientific_names, on="tax_id", how="left")

# 提取 family 和 genus
families = taxonomy[taxonomy["rank"] == "family"][["tax_id", "name_txt"]].rename(columns={"name_txt": "family"})
genera = taxonomy[taxonomy["rank"] == "genus"][["tax_id", "parent_tax_id", "name_txt"]].rename(columns={"name_txt": "genus"})

# 建立 tax_id -> parent_tax_id, rank 映射
parent_map = taxonomy.set_index("tax_id")[["parent_tax_id", "rank"]].to_dict(orient="index")

def find_family(tax_id):
    """递归查找给定 tax_id 的上级 family"""
    visited = set()
    while tax_id in parent_map and tax_id not in visited:
        visited.add(tax_id)
        info = parent_map[tax_id]
        if info["rank"] == "family":
            return tax_id
        tax_id = info["parent_tax_id"]
    return None

# 给 genus 找到对应 family
genera["family_tax_id"] = genera["parent_tax_id"].apply(find_family)

# 过滤掉没有找到 family 的 genus
genera_family = genera.dropna(subset=["family_tax_id"])

# 合并 family 名称
genera_family = genera_family.merge(
    families,
    left_on="family_tax_id",
    right_on="tax_id",
    how="left"
)[["family", "genus"]]

# 分组统计
family_grouped_count = genera_family.groupby("family")["genus"].agg(
    genera_list=lambda x: ", ".join(sorted(set(x))),
    genera_count=lambda x: len(set(x))
).reset_index()

# 导出 Excel
output_file = "~/Desktop/family_genus_list.xlsx"
family_grouped_count.to_excel(output_file, index=False)

print(f"结果已保存为 {output_file}")

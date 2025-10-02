# There is likely a better solution to this, but for now this is how content will be stored
def outputContent(topic):
    if topic == "GT":
        return {"Graph" : "Definition for this object", "Neighbor Set":"Definition for this object", "Degree":"Definition for this object",
                "Complete Graph":"Definition for this object"}
    elif topic == "LA":
        return {"Matrix": "stuff", "Vector space": "other stuff", "Orthogonality":"some definition involving inner products and other stuff like that"}
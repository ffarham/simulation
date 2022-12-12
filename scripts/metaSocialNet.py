import facebook
import json

# custom python package containing all API tokens
from tokens import tokens

def main():
    graph = facebook.GraphAPI(tokens.getMeta())
    profile = graph.get_object("me")

    print(json.dumps(profile, indent=4))

if __name__ == '__main__':
    main()
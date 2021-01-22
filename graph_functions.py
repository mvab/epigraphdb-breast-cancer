def start_graph_session(env):
    
    from neo4j import GraphDatabase
    env.read_env()

    EPIGRAPHDB_SERVER = env.str("EPIGRAPHDB_SERVER")
    EPIGRAPHDB_USER = "neo4j"
    EPIGRAPHDB_PORT = env.str("EPIGRAPHDB_PORT")
    EPIGRAPHDB_PASSWORD = env.str("EPIGRAPHDB_PASSWORD")
    print(EPIGRAPHDB_PORT)

    epigraphdb_driver = GraphDatabase.driver(
        "bolt://{server_name}:{port}".format(
            server_name=EPIGRAPHDB_SERVER, port=EPIGRAPHDB_PORT),
        auth=(EPIGRAPHDB_USER, EPIGRAPHDB_PASSWORD))
    
    session = epigraphdb_driver.session()
    return session



def query_to_df(session, query, print_query=True):
    import pandas as pd
    if print_query:
        print(query)
    data=session.run(query).data()
    df = pd.json_normalize(data)
    return(df)